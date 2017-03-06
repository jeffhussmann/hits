import bokeh
import numpy as np
import glob
import os.path
import pandas as pd
from itertools import cycle
from . import external_coffeescript, colors_list
from meta import ToggleLegend
from bokeh.models.annotations import LegendItem

def load_data(exp_dirs):
    def read_csv(fn):
        fh = open(fn)
        line = ''
        while not line.startswith('Sample Name'):
            line = fh.readline()
            
        name = line.strip().split(',', 1)[1].translate(None, '"')
        
        while not line.startswith('Time,Value'):
            line = fh.readline()
        
        times = []
        values = []
        
        line = fh.readline().strip()
        while line:
            time, value = map(float, line.strip().split(','))
            times.append(time)
            values.append(value)
            
            line = fh.readline().strip()
            
        series = pd.Series(values, index=times, name=name)
        
        return series

    ordered = ['Sample{0}'.format(i) for i in range(1, 12)] + ['Ladder']
    name_to_order = {n: i + 1 for i, n in enumerate(ordered)}
    def fn_to_order(fn):
        root, ext = os.path.splitext(fn)
        name = root.split('_')[-1]
        if name == 'Results':
            return None

        order = name_to_order[name]
        return order

    exps = {
        'raw': {},
        'nonmarker area': {},
        'marker area': {},
        'descriptions': {},
    }

    for exp_dir in exp_dirs:
        try:
            description_fn = '{0}/descriptions.txt'.format(exp_dir)
            descriptions = dict(line.strip().split('\t') for line in open(description_fn))
            descriptions['Ladder'] = 'ladder'
        except IOError:
            descriptions = {}

        head, tail = os.path.split(exp_dir)
        for key in exps:
            exps[key][tail] = {}
        fns = glob.glob('{0}/*.csv'.format(exp_dir))
        for fn in fns:
            data = read_csv(fn)
            order = fn_to_order(fn)
            if order is None:
                continue
            name = '{0}:{1:02d}_{2}'.format(tail, order, data.name)

            marker_area = sum(data[22:23])
            nonmarker_area = sum(data[23:])

            exps['raw'][tail][name] = data
            exps['nonmarker area'][tail][name] = data / nonmarker_area
            exps['marker area'][tail][name] = data / marker_area
            exps['descriptions'][tail][name] = descriptions.get(data.name, data.name)
        
    return exps

def extract_ladder_peak_times(series):
    ladder_peaks = [
        ((20, 25), 25),
        ((25, 30), 200),
        ((30, 35), 500),
        ((35, 40), 1000),
        ((40, 46), 2000),
        ((46, 52.5), 4000),
    ]
    peak_times = {series[start:end].idxmax(): nt for (start, end), nt in ladder_peaks}
    return peak_times

def plot(exps,
         initial_sub_group_selections=None,
         unselected_alpha=0.2,
        ):
    if initial_sub_group_selections is None:
        initial_sub_group_selections = []
    groupings = {group_name: sorted(exps['raw'][group_name]) for group_name in exps['raw']}

    sources = {}
    colors = {}
    color_iter = cycle(colors_list)

    highest_level_keys = [k for k in sorted(exps) if k != 'descriptions']

    for key in ['plotted'] + highest_level_keys:
        if key == 'plotted':
            normalization = 'raw'
        else:
            normalization = key
        
        sources[key] = {}
        for group_name in sorted(exps[normalization]):
            for sample_name in sorted(exps[normalization][group_name]):
                if sample_name not in colors:
                    colors[sample_name] = color_iter.next()
                series = exps[normalization][group_name][sample_name]
                source = bokeh.models.ColumnDataSource()
                source.data['x'] = list(series.index)
                source.data['y'] = list(series)
                source.data['name'] = [sample_name] * len(series)
                source.name = 'source_{0}_{1}'.format(sample_name, key)
                sources[key][sample_name] = source
    
    sample_names = []
    descriptions = []
    for group_name, group in sorted(exps['descriptions'].items()):
        for sample_name, description in sorted(group.items()):
            sample_names.append(sample_name)
            descriptions.append(description)

    full_source = bokeh.models.ColumnDataSource(name='full_source')
    full_source.data = {
        'sample_name': sample_names,
        'description': descriptions,
    }

    tools = [
        'pan',
        'tap',
        'box_zoom',
        'wheel_zoom',
        'save',
        'reset',
        'undo',
    ]

    x_range = bokeh.models.Range1d(17, 70, bounds=(17, 70))
    fig = bokeh.plotting.figure(plot_width=1200,
                                plot_height=600,
                                tools=tools, active_scroll='wheel_zoom',
                                x_range=x_range,
                               )

    random_group = exps['raw'].keys()[0]
    series = exps['raw'][random_group]['{0}:12_Ladder'.format(random_group)]
    peak_times = extract_ladder_peak_times(series)
    
    for time, nt in peak_times.items():
        line = bokeh.models.annotations.Span(location=time,
                                             dimension='height',
                                             line_color='black',
                                             line_alpha=0.8,
                                             line_dash='dashed',
                                            )
        fig.renderers.append(line)

    fig.xgrid.grid_line_color = None
    fig.ygrid.grid_line_color = None

    convert_tick = '''
    peak_times = {dict};
    return peak_times[tick];
    '''.format(dict=peak_times)
    fig.xaxis.formatter = bokeh.models.FuncTickFormatter(code=convert_tick)
    fig.xaxis.ticker = bokeh.models.FixedTicker(ticks=list(peak_times))

    all_legend_items = []
    initial_legend_items = []
    lines = []
    for sample_name, source in sources['plotted'].items():
        if sample_name in initial_sub_group_selections:
            color = colors[sample_name]
            line_width = 2
            line_alpha = 0.95
        else:
            color = 'black'
            line_width = 1
            if len(initial_sub_group_selections) > 0:
                line_alpha = unselected_alpha
            else:
                line_alpha = 0.6

        line = fig.line(x='x',
                        y='y',
                        color=color,
                        source=source,
                        line_width=line_width,
                        line_alpha=line_alpha,
                        line_join='round',
                        nonselection_line_color=colors[sample_name],
                        nonselection_line_alpha=unselected_alpha,
                        hover_alpha=1.0,
                        hover_color=colors[sample_name],
                       )

        line.hover_glyph.line_width = 4
        line.name = 'line_{0}'.format(sample_name)
        lines.append(line)

        legend_item = LegendItem(label=sample_name, renderers=[line])
        all_legend_items.append(legend_item)
        if sample_name in initial_sub_group_selections:
            initial_legend_items.append(legend_item)
        
    legend = ToggleLegend(name='legend',
                          items=initial_legend_items,
                          all_items=all_legend_items,
                         )
    fig.add_layout(legend)
    
    source_callback = external_coffeescript('bioanalyzer_selection',
                                            args=dict(full_source=full_source),
                                           )
    for source in sources['plotted'].values():
        source.callback = source_callback

    hover = bokeh.models.HoverTool(line_policy='interp',
                                   renderers=lines,
                                  )
    hover.tooltips = [('name', '@name')]
    fig.add_tools(hover)

    sub_group_callback = external_coffeescript('bioanalyzer_sub_group',
                                               args=dict(full_source=full_source),
                                              )

    top_group_callback = external_coffeescript('bioanalyzer_top_group',
                                               args=dict(full_source=full_source),
                                              )

    top_groups = []
    sub_groups = []
    for top_name, sub_names in sorted(groupings.items()):
        width = 75 + max(len(l) for l in sub_names) * 6
        top = bokeh.models.widgets.CheckboxGroup(labels=[top_name],
                                                 active=[],
                                                 width=width,
                                                 name='top_{0}'.format(top_name),
                                                 callback=top_group_callback,
                                                )
        top_groups.append(top)

        sub_active = [i for i, n in enumerate(sub_names) if n in initial_sub_group_selections]
        sub = bokeh.models.widgets.CheckboxGroup(labels=sub_names,
                                                 active=sub_active,
                                                 width=width,
                                                 callback=sub_group_callback,
                                                 name='sub_{0}'.format(top_name),
                                                )
        sub_groups.append(sub)

    
    clear_selection = bokeh.models.widgets.Button(label='Clear selection')
    clear_selection.callback = external_coffeescript('bioanalyzer_clear_selection',
                                                     args=dict(full_source=full_source),
                                                    )
    
    highest_level_chooser = bokeh.models.widgets.Select(options=highest_level_keys,
                                                        value=highest_level_keys[2],
                                                       )
    callback_name = 'metacodon_highest_level'
    
    injection_sources = []
    for key in highest_level_keys:
        injection_sources.extend(sources[key].values())
    injection = {'ensure_no_collision_{0}'.format(i): v for i, v in enumerate(injection_sources)}

    highest_level_chooser.callback = external_coffeescript(callback_name,
                                                           args=injection,
                                                          )
    table_col_names = [
        ('sample_name', 250),
        ('description', 850),
    ]

    columns = []
    for col_name, width in table_col_names:
        column = bokeh.models.widgets.TableColumn(field=col_name,
                                                  title=col_name,
                                                  formatter=None,
                                                  width=width,
                                                 )
        columns.append(column)

    filtered_data = {k: [] for k, _ in table_col_names}
    filtered_source = bokeh.models.ColumnDataSource(data=filtered_data, name='table_source')

    table = bokeh.models.widgets.DataTable(source=filtered_source,
                                           columns=columns,
                                           width=1200,
                                           height=400,
                                           sortable=False,
                                           name='table',
                                           row_headers=False,
                                          )
    grid = bokeh.layouts.layout([
        sub_groups + [bokeh.layouts.widgetbox([highest_level_chooser, clear_selection])],
        [fig],
        [table],
    ])

    bokeh.io.show(grid)
