import bokeh
from bokeh.core.properties import List, Instance
from bokeh.models.tickers import SingleIntervalTicker
from bokeh.models.annotations import LegendItem
import numpy as np
import positions
from itertools import cycle
from . import external_coffeescript, colors_list

class ToggleLegend(bokeh.models.annotations.Legend):
    all_items = List(Instance(LegendItem))

    __implementation__ = '''
import {Legend} from "models/annotations/legend"
import * as p from "core/properties"

export class ToggleLegend extends Legend
    type: "ToggleLegend"
    @define {
        all_items: [p.Array, []]
    }
'''

def codon(xs, ys, colors, groupings,
          toggle_resolution=True,
          y_max=5,
          x_lims=(-25, 25),
          unselected_alpha=0.2,
          initial_menu_selection=None,
          initial_top_group_selections=None,
          initial_sub_group_selections=None,
          intial_resolution='codon',
         ):
    if initial_top_group_selections is None:
        initial_top_group_selections = []
    if initial_sub_group_selections is None:
        initial_sub_group_selections = []

    initial_sub_group_selections = set(initial_sub_group_selections)
    for top_name, sub_names in sorted(groupings.items()):
        if top_name in initial_top_group_selections:
            initial_sub_group_selections.update(sub_names)
    
    sources = {}

    highest_level_keys = sorted(ys)
    
    menu_options = sorted(ys[highest_level_keys[0]].values()[0].keys())

    if initial_menu_selection is None:
        initial_menu_selection = menu_options[0]
    if initial_menu_selection not in menu_options:
        raise ValueError('{0} not in {1}'.format(initial_menu_selection, menu_options))

    for key in ['plotted'] + highest_level_keys:
        if key == 'plotted':
            resolution = intial_resolution
        else:
            resolution = key
        
        sources[key] = {}
        for checkbox_name in sorted(ys[resolution]):
            source = bokeh.models.ColumnDataSource(ys[resolution][checkbox_name])
            source.data['x'] = xs[resolution]
            source.data['y'] = source.data[initial_menu_selection]
            source.data['name'] = [checkbox_name] * len(xs[resolution])
            source.name = 'source_{0}_{1}'.format(checkbox_name, key)
            sources[key][checkbox_name] = source
   
    tools = [
        'pan',
        'tap',
        'box_zoom',
        'wheel_zoom',
        'save',
        'reset',
        'undo',
    ]

    fig = bokeh.plotting.figure(plot_width=1200,
                                plot_height=800,
                                tools=tools,
                                active_scroll='wheel_zoom',
                                name='figure',
                               )

    fig.grid.grid_line_alpha = 0.4

    fig.y_range = bokeh.models.Range1d(0, y_max)
    fig.y_range.name = 'y_range'
    fig.x_range = bokeh.models.Range1d(*x_lims)
    fig.x_range.name = 'x_range'

    range_callback = external_coffeescript('metacodon_range')
    fig.y_range.callback = range_callback
    fig.x_range.callback = range_callback

    fig.xaxis.axis_label = 'Offset ({0}s)'.format(intial_resolution)
    fig.yaxis.axis_label = 'Mean relative enrichment'

    fig.xaxis.name = 'x_axis'
    fig.yaxis.name = 'y_axis'

    fig.xaxis.axis_label_text_font_style = 'normal'
    fig.yaxis.axis_label_text_font_style = 'normal'

    legend_items = []
    initial_legend_items = []
    lines = []
    for checkbox_name, source in sources['plotted'].items():
        if checkbox_name in initial_sub_group_selections:
            color = colors[checkbox_name]
            line_width = 2
            line_alpha = 0.95
            circle_visible = True
        else:
            color ='black'
            line_width = 1
            circle_visible = False
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
                        nonselection_line_color=colors[checkbox_name],
                        nonselection_line_alpha=unselected_alpha,
                        hover_alpha=1.0,
                        hover_color=colors[checkbox_name],
                       )
        line.hover_glyph.line_width = 3
        line.name = 'line_{0}'.format(checkbox_name)
        lines.append(line)
        
        circle = fig.circle(x='x',
                            y='y',
                            color=colors[checkbox_name],
                            source=source,
                            size=2.5,
                            fill_alpha=0.95,
                            line_alpha=0.95,
                            visible=circle_visible,
                            hover_alpha=0.95,
                            hover_color=colors[checkbox_name],
                           )
        circle.hover_glyph.visible = True
        circle.name = 'circle_{0}'.format(checkbox_name)
    
        legend_item = LegendItem(label=checkbox_name, renderers=[line])
        legend_items.append(legend_item)
        if checkbox_name in initial_sub_group_selections:
            initial_legend_items.append(legend_item)
        
    legend = ToggleLegend(name='legend',
                          items=initial_legend_items,
                          all_items=legend_items,
                         )
    fig.add_layout(legend)
    fig.legend.location = 'top_right'
    fig.legend.background_fill_alpha = 0.5

    fig.xaxis[0].ticker = bokeh.models.tickers.SingleIntervalTicker(interval=3, num_minor_ticks=3)

    source_callback = external_coffeescript('metacodon_selection')
    for source in sources['plotted'].values():
        source.callback = source_callback

    hover = bokeh.models.HoverTool(line_policy='interp',
                                   renderers=lines,
                                  )
    hover.tooltips = [('name', '@name')]
    fig.add_tools(hover)

    zero_x = bokeh.models.annotations.Span(location=0,
                                           dimension='height',
                                           line_color='black',
                                           line_alpha=0.8,
                                          )
    fig.renderers.append(zero_x)

    one_y = bokeh.models.annotations.Span(location=1,
                                           dimension='width',
                                           line_color='black',
                                           line_alpha=0.8,
                                          )
    fig.renderers.append(one_y)

    menu = bokeh.models.widgets.MultiSelect(options=menu_options,
                                            value=[initial_menu_selection],
                                            size=min(30, len(menu_options)),
                                           )
    menu.callback = external_coffeescript('metacodon_menu')

    sub_group_callback = external_coffeescript('metacodon_sub_group',
                                               format_kwargs=dict(color_unselected='false'),
                                              )

    top_group_callback = external_coffeescript('metacodon_top_group')

    top_groups = []
    sub_groups = []

    for top_name, sub_names in sorted(groupings.items()):
        
        width = 75 + max(len(l) for l in sub_names) * 6

        top_active = [0] if top_name in initial_top_group_selections else []
        top = bokeh.models.widgets.CheckboxGroup(labels=[top_name],
                                                 active=top_active,
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

    if toggle_resolution:
        highest_level_chooser = bokeh.models.widgets.RadioGroup(labels=['codon resolution', 'nucleotide resolution'],
                                                                active=0 if intial_resolution == 'codon' else 1,
                                                                name='highest_level_chooser',
                                                               )
        callback_name = 'metacodon_resolution'
    else:
        highest_level_chooser = bokeh.models.widgets.Select(options=highest_level_keys,
                                                            value=highest_level_keys[0],
                                                            name='highest_level_chooser',
                                                           )
        callback_name = 'metacodon_highest_level'
    
    injection_sources = []
    for key in highest_level_keys:
        injection_sources.extend(sources[key].values())
    injection = {'ensure_no_collision_{0}'.format(i): v for i, v in enumerate(injection_sources)}

    highest_level_chooser.callback = external_coffeescript(callback_name,
                                                           args=injection,
                                                          )

    clear_selection = bokeh.models.widgets.Button(label='Clear selection')
    clear_selection.callback = external_coffeescript('metacodon_clear_selection')

    grid = [
        top_groups,
        sub_groups,
        [fig, bokeh.layouts.widgetbox([menu, highest_level_chooser, clear_selection])],
    ]

    bokeh.io.show(bokeh.layouts.layout(grid))

def gene(densities,
         groupings,
         y_max=5,
         unselected_alpha=0.2,
         assignment='three_prime',
         min_density=0,
         start_grey=False,
         initial_resolution='nucleotide',
         initial_top_group_selections=None,
         initial_sub_group_selections=None,
        ):
    if initial_top_group_selections is None:
        initial_top_group_selections = []
    if initial_sub_group_selections is None:
        initial_sub_group_selections = []

    initial_sub_group_selections = set(initial_sub_group_selections)
    for top_name, sub_names in sorted(groupings.items()):
        if top_name in initial_top_group_selections:
            initial_sub_group_selections.update(sub_names)
    
    exp_names = sorted(densities['codon'])
    sources = {}

    highest_level_keys = sorted(densities)
    
    max_before = 90
    max_after = 250
    xs_dict = {
        'codon': {
            'start_codon': np.arange(-max_before, max_after),
            'stop_codon': np.arange(-max_after, max_before),
        },
        'nucleotide': {
            'start_codon': np.arange(-3 * max_before, 3 * max_after),
            'stop_codon': np.arange(-3 * max_after, 3 * max_before),
        },
    }

    for key in ['plotted'] + highest_level_keys:
        if key == 'plotted':
            resolution = initial_resolution
        else:
            resolution = key
        
        sources[key] = {}
        for exp_name in exp_names:
            source = bokeh.models.ColumnDataSource()
            for landmark in ['start_codon', 'stop_codon']:
                xs = xs_dict[resolution][landmark]

                if resolution == 'codon':
                    actual_assignment = assignment
                else:
                    actual_assignment = 'three_prime'

                density_type = positions.MetageneDensityType(actual_assignment,
                                                             landmark,
                                                             'all',
                                                             'none_nearby',
                                                             min_density,
                                                            )
                ys = densities[resolution][exp_name][str(density_type)]['means'][landmark, xs]

                source.data['xs_{0}'.format(landmark)] = xs
                source.data['ys_{0}'.format(landmark)] = ys
                source.data['name'] = [exp_name] * len(xs)
            
            source.name = 'source_{0}_{1}'.format(exp_name, key)
            sources[key][exp_name] = source

    initial_before = 20
    initial_after = 150

    initial_lims = {
        'start_codon': (-initial_before, initial_after),
        'stop_codon': (-initial_after, initial_before),
    }

    x_ranges = {}
    x_range_callback = external_coffeescript('metagene_x_range')
    for landmark in ('start_codon', 'stop_codon'):
        x_range = bokeh.models.Range1d(*initial_lims[landmark])
        x_range.name = 'x_range_{0}'.format(landmark)
        x_range.callback = x_range_callback
        x_ranges[landmark] = x_range

    y_range = bokeh.models.Range1d(0, y_max)
    y_range.name = 'y_range'
    y_range.callback = external_coffeescript('metagene_y_range')

    tools = [
        'pan',
        'tap',
        'box_zoom',
        'wheel_zoom',
        'save',
        'reset',
        'undo',
    ]

    figs = {}
    for key, y_axis_location in [('start_codon', 'left'), ('stop_codon', 'right')]:
        figs[key] = bokeh.plotting.figure(plot_width=800,
                                          plot_height=600,
                                          x_range=x_ranges[key],
                                          y_range=y_range,
                                          y_axis_location=y_axis_location,
                                          tools=tools,
                                          active_scroll='wheel_zoom',
                                         )

    lines = {
        'start_codon': [],
        'stop_codon': [],
    }

    legend_items = []
    initial_legend_items = []

    colors = dict(zip(exp_names, cycle(colors_list)))
    for exp_name in exp_names:
        if exp_name in initial_sub_group_selections:
            color = colors[exp_name]
            line_width = 2
            line_alpha = 0.95
            circle_visible = True
        else:
            color = 'black'
            line_width = 1
            circle_visible = False
            if len(initial_sub_group_selections) > 0:
                line_alpha = unselected_alpha
            else:
                line_alpha = 0.6

        for landmark in ('start_codon', 'stop_codon'):
            line = figs[landmark].line(x='xs_{0}'.format(landmark),
                                       y='ys_{0}'.format(landmark),
                                       source=sources['plotted'][exp_name],
                                       color=color,
                                       nonselection_line_color=colors[exp_name],
                                       nonselection_line_alpha=unselected_alpha,
                                       hover_alpha=1.0,
                                       hover_color=colors[exp_name],
                                       line_width=line_width,
                                       line_alpha=line_alpha,
                                      )
            line.hover_glyph.line_width = 4
            line.name = 'line_{0}'.format(exp_name)
            lines[landmark].append(line)
            
            circle = figs[landmark].circle(x='xs_{0}'.format(landmark),
                                           y='ys_{0}'.format(landmark),
                                           source=sources['plotted'][exp_name],
                                           size=2,
                                           color=colors[exp_name],
                                           fill_alpha=0.9,
                                           line_alpha=0.9,
                                           visible=circle_visible,
                                           hover_alpha=1.0,
                                           hover_color=colors[exp_name],
                                          )
            circle.hover_glyph.visible = True
            circle.name = 'circle_{0}'.format(exp_name)

            if landmark == 'stop_codon':
                legend_item = LegendItem(label=exp_name, renderers=[line])
                legend_items.append(legend_item)
                if exp_name in initial_sub_group_selections:
                    initial_legend_items.append(legend_item)

    legend = ToggleLegend(name='legend',
                          items=initial_legend_items,
                          all_items=legend_items,
                         )
    figs['stop_codon'].add_layout(legend)
    figs['stop_codon'].legend.location = 'top_left'
    figs['stop_codon'].legend.background_fill_alpha = 0.5
    
    for landmark, fig in figs.items():
        zero_x = bokeh.models.annotations.Span(location=0,
                                               dimension='height',
                                               line_color='black',
                                               line_alpha=0.5,
                                              )
        fig.renderers.append(zero_x)

        one_y = bokeh.models.annotations.Span(location=1,
                                               dimension='width',
                                               line_color='black',
                                               line_alpha=0.5,
                                              )
        fig.renderers.append(one_y)
    
    source_callback = external_coffeescript('metacodon_selection')
    for source in sources['plotted'].values():
        source.callback = source_callback

    for landmark in ['start_codon', 'stop_codon']:
        hover = bokeh.models.HoverTool(line_policy='interp',
                                       renderers=lines[landmark],
                                      )
        hover.tooltips = [('name', '@name')]
        figs[landmark].add_tools(hover)

    resolution = bokeh.models.widgets.RadioGroup(labels=['codon resolution', 'nucleotide resolution'],
                                                 active=0 if initial_resolution == 'codon' else 1,
                                                )
    resolution.name = 'resolution'
    
    injection_sources = []
    for key in highest_level_keys:
        injection_sources.extend(sources[key].values())
    injection = {'ensure_no_collision_{0}'.format(i): v for i, v in enumerate(injection_sources)}

    resolution.callback = external_coffeescript('metacodon_resolution',
                                                args=injection,
                                               )
    
    sub_group_callback = external_coffeescript('metacodon_sub_group',
                                               format_kwargs=dict(color_unselected='false'),
                                              )

    top_group_callback = external_coffeescript('metacodon_top_group')

    top_groups = []
    sub_groups = []
    width = 100
    for top_name, sub_names in sorted(groupings.items()):
        width = 75 + max(len(l) for l in sub_names) * 6
        
        top_active = [0] if top_name in initial_top_group_selections else []
        top = bokeh.models.widgets.CheckboxGroup(labels=[top_name],
                                                 active=top_active,
                                                 width=width,
                                                 name='top_{0}'.format(top_name),
                                                 callback=top_group_callback,
                                                )
        top_groups.append(top)

        sub_active = [i for i, n in enumerate(sub_names) if n in initial_sub_group_selections]
        sub = bokeh.models.widgets.CheckboxGroup(labels=sorted(sub_names),
                                                 active=sub_active,
                                                 width=width,
                                                 callback=sub_group_callback,
                                                 name='sub_{0}'.format(top_name),
                                                )
        sub_groups.append(sub)

    plots = bokeh.layouts.gridplot([[figs['start_codon'], figs['stop_codon']]])
    grid = [
        top_groups,
        sub_groups,
        [plots, resolution],
    ]
    bokeh.io.show(bokeh.layouts.layout(grid))

def lengths(ys, group_by='experiment', groupings=None,
            initial_top_group_selections=[],
            initial_sub_group_selections=[],
           ):
    sources = {}

    first_exp_name = ys.keys()[0]

    # ys has experiments as top level keys and types as lower level keys.
    # Need a dictionary with checkbox names as the top keys, so if grouping by
    # type, don't need to do anything. If grouping by experiment, need to invert
    # the order.

    if group_by == 'experiment':
        menu_options = sorted(ys)
        checkbox_names = sorted(ys[first_exp_name])
        
        rekeyed_ys = {}
        for checkbox_name in checkbox_names:
            rekeyed_ys[checkbox_name] = {}
            for menu_name in menu_options:
                rekeyed_ys[checkbox_name][menu_name] = ys[menu_name][checkbox_name]

        ys = rekeyed_ys
    elif group_by == 'type':
        menu_options = sorted(ys[first_exp_name])
        checkbox_names = sorted(ys)

    for checkbox_name in ys:
        for menu_name in ys[checkbox_name]:
            ys[checkbox_name][menu_name] = ys[checkbox_name][menu_name][:200]
    
    if groupings is None:
        groupings = {n: [n] for n in checkbox_names}

    colors = dict(zip(checkbox_names, cycle(colors_list)))

    initial_menu_selection = None
    if initial_menu_selection is None:
        initial_menu_selection = menu_options[0]

    if initial_menu_selection not in menu_options:
        raise ValueError('{0} not in {1}'.format(initial_menu_selection, menu_options))

    for checkbox_name in checkbox_names:
        source = bokeh.models.ColumnDataSource(ys[checkbox_name])

        xs = range(len(ys[checkbox_name][menu_options[0]]))
        source.data['x'] = xs

        source.data['y'] = source.data[initial_menu_selection]

        source.data['name'] = [checkbox_name] * len(xs)
        source.name = 'source_{0}_plotted'.format(checkbox_name)
        sources[checkbox_name] = source

    tools = [
        'pan',
        'tap',
        'box_zoom',
        'wheel_zoom',
        'save',
        'reset',
        'undo',
    ]

    fig = bokeh.plotting.figure(plot_width=1200,
                                plot_height=800,
                                tools=tools,
                                active_scroll='wheel_zoom',
                                name='figure',
                               )

    fig.grid.grid_line_alpha = 0.4

    #fig.y_range = bokeh.models.Range1d(0, y_max)
    fig.y_range.name = 'y_range'
    #fig.x_range = bokeh.models.Range1d(-x_max, x_max)
    fig.x_range.name = 'x_range'

    #range_callback = external_coffeescript('metacodon_range')
    #fig.y_range.callback = range_callback
    #fig.x_range.callback = range_callback

    fig.xaxis.axis_label = 'Length'
    fig.yaxis.axis_label = 'y'

    fig.xaxis.name = 'x_axis'
    fig.yaxis.name = 'y_axis'

    fig.xaxis.axis_label_text_font_style = 'normal'
    fig.yaxis.axis_label_text_font_style = 'normal'

    legend_items = []
    initial_legend_items = []
    lines = []
    for checkbox_name, source in sources.items():
        if checkbox_name in initial_sub_group_selections:
            color = colors[checkbox_name]
            line_width = 2
            line_alpha = 0.95
            circle_visible = True
        else:
            color = colors[checkbox_name]
            line_width = 1
            circle_visible = False
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
                        nonselection_line_color=colors[checkbox_name],
                        nonselection_line_alpha=0.5,
                        hover_alpha=1.0,
                        hover_color=colors[checkbox_name],
                       )
        line.hover_glyph.line_width = 3
        line.name = 'line_{0}'.format(checkbox_name)
        lines.append(line)
        
        circle = fig.circle(x='x',
                            y='y',
                            color=colors[checkbox_name],
                            source=source,
                            size=2,
                            fill_alpha=0.95,
                            line_alpha=0.95,
                            visible=circle_visible,
                            hover_alpha=0.95,
                            hover_color=colors[checkbox_name],
                           )
        circle.hover_glyph.visible = True
        circle.name = 'circle_{0}'.format(checkbox_name)
    
        legend_item = LegendItem(label=checkbox_name, renderers=[line])
        legend_items.append(legend_item)
        if checkbox_name in initial_sub_group_selections:
            initial_legend_items.append(legend_item)
        
    legend = ToggleLegend(name='legend',
                          items=initial_legend_items,
                          all_items=legend_items,
                         )
    fig.add_layout(legend)
    
    source_callback = external_coffeescript('metacodon_selection')
    for source in sources.values():
        source.callback = source_callback
    
    hover = bokeh.models.HoverTool(line_policy='interp',
                                   renderers=lines,
                                  )
    hover.tooltips = [('name', '@name')]
    fig.add_tools(hover)
    
    menu = bokeh.models.widgets.MultiSelect(options=menu_options,
                                            value=[initial_menu_selection],
                                            size=min(40, len(menu_options)),
                                           )
    menu.callback = external_coffeescript('metacodon_menu')
    
    sub_group_callback = external_coffeescript('metacodon_sub_group',
                                               format_kwargs=dict(color_unselected='true'),
                                              )
    top_group_callback = external_coffeescript('metacodon_top_group')

    top_groups = []
    sub_groups = []
    width = 100
    for top_name, sub_names in sorted(groupings.items()):
        width = 75 + max(len(l) for l in sub_names) * 6
        
        top_active = [0] if top_name in initial_top_group_selections else []
        top = bokeh.models.widgets.CheckboxGroup(labels=[top_name],
                                                 active=top_active,
                                                 width=width,
                                                 name='top_{0}'.format(top_name),
                                                 callback=top_group_callback,
                                                )
        top_groups.append(top)

        sub_active = [i for i, n in enumerate(sub_names) if n in initial_sub_group_selections]
        sub = bokeh.models.widgets.CheckboxGroup(labels=sorted(sub_names),
                                                 active=sub_active,
                                                 width=width,
                                                 callback=sub_group_callback,
                                                 name='sub_{0}'.format(top_name),
                                                )
        sub_groups.append(sub)

    grid = [
        top_groups,
        sub_groups,
        [fig, menu],
    ]
    bokeh.io.show(bokeh.layouts.layout(grid))
