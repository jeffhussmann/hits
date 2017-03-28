import numpy as np
import bokeh
import bokeh.io
import bokeh.plotting
import pandas as pd
import matplotlib.colors
import matplotlib.cm
import os.path
import glob
import IPython.display
from collections import defaultdict
from external_coffeescript import build_callback
import meta

bokeh.io.output_notebook()

def build_selected(indices):
    pvd = bokeh.core.property.containers.PropertyValueDict
    pvl = bokeh.core.property.containers.PropertyValueList

    selected = pvd({
        '0d': pvd({
            'glyph': None,
            'indices': pvl(),
        }),
        '1d': pvd({
            'indices': pvl(indices),
        }),
        '2d': pvd(),
    })

    return selected

def scatter(df=None,
            hover_keys=None,
            table_keys=None,
            color_by=None,
            label_by=None,
            size=800,
            axis_label_size=20,
            log_scale=False,
            volcano=False,
            heatmap=False,
            grid=False,
            marker_size=6,
            initial_selection=None,
            initial_xy_names=None,
            data_lims=None,
            alpha_widget_type='slider',
            hide_widgets=None,
            identical_bins=True,
            num_bins=100,
           ):
    ''' Makes an interactive scatter plot using bokeh. Call without any
    arguments for an example using data from Jan et al. Science 2014.

    Args:
            df: A pandas DataFrame with columns containing numerical data to
                plot. 
                Index values will be used as labels for points.
                Any text columns will be searchable through the 'Search:' field.
                Any boolean columns will be used to define subsets of points for
                selection from a dropdown menu.
                If df is None, loads example data from Jan et al. Science 2014.

            hover_keys: Names of columns in df to display in the tooltip that
                appears when you hover over a point.

            table_keys: Names of columns in df to display in the table below the
                plot that is populated with the selected points from the figure.

            color_by: The name of a column in df to use as colors of points, or
                a list of such names to choose from a menu.
            
            label_by: The name of a column in df to use as labels of points, or
                a list of such names to choose from a menu. If None, df.index
                is used.

            size: Size of the plot in pixels.

            marker_size: Size of the scatter circles.

            heatmap: If True, displays a heatmap of correlations between
                numerical columns in df that can be clicked to select columns
                to scatter.

            grid: If True, defaults to grid instead of diagonal landmarks.

            volcano: If True, make some tweaks suitable for volcano plots.

            log_scale: If not False, plot on a log scale with base 10 (or, if a
                set to a number, with base log_scale.)

            axis_label_size: Size of the font used for axis labels.

            initial_selection: Names of index value to initially highlight.

            initial_xy_names: Tuple (x_name, y_name) of datasets to initially
                display on x- and y-axes.

            alpha_widget_type: Type of widget ('slider' or 'text') to control
                scatter circle transparency.

            hide_widgets: List of widgets to not display. Possible options are
                ['table', 'alpha', 'marker_size', 'search', 'subset_menu',
                 'grid_radio_buttons'].

            identical_bins: If True, use the same set of bins for histograms of
                all data sets. If False, use bins tailored to the range of each
                data set.

            num_bins: Number of bins to use for marginal histograms.
    '''

    if hover_keys is None:
        hover_keys = []

    if table_keys is None:
        table_keys = []
    
    if hide_widgets is None:
        hide_widgets = []

    if volcano:
        grid = True

    if df is None:
        # Load example data.
        fn = os.path.join(os.path.dirname(__file__), 'example_df.txt')
        df = pd.read_csv(fn, index_col='alias')

        # Override some arguments.
        log_scale = True
        hover_keys = ['systematic_name', 'short_description']
        table_keys = ['systematic_name', 'description']
        grid = False
        heatmap = True

    # Infer column types.
    
    scatter_source = bokeh.models.ColumnDataSource(data=df, name='scatter_source')

    if 'index' in scatter_source.data:
        scatter_source.data['_index'] = scatter_source.data['index']

    if df.index.name is None:
        df.index.name = 'index'

    if initial_selection is None:
        initial_selection = []

    initial_indices = [i for i, n in enumerate(df.index) if n in initial_selection]

    numerical_cols = [n for n in df.columns if df[n].dtype in [float, int]]

    object_cols = [n for n in df.columns if df[n].dtype is np.dtype('O')]
    if df.index.dtype is np.dtype('O'):
        object_cols.append(df.index.name)

    bool_cols = [n for n in df.columns if df[n].dtype is np.dtype('bool')]

    subset_indices = {n: [i for i, v in enumerate(df[n]) if v] for n in bool_cols}

    # Set up the actual scatter plot.
    
    tools = [
        'reset',
        'undo',
        'pan',
        'box_zoom',
        'box_select',
        'tap',
        'wheel_zoom',
        'save',
    ]
    
    fig_kwargs = dict(
        plot_width=size,
        plot_height=size,
        tools=tools,
        lod_threshold=10000,
        name='scatter_fig',
    )

    min_border = 80

    if log_scale:
        if log_scale is True:
            log_scale = 10
        fig_kwargs['y_axis_type'] = 'log'
        fig_kwargs['x_axis_type'] = 'log'
    
    fig = bokeh.plotting.figure(**fig_kwargs)
    fig.toolbar.logo = None
    fig.toolbar_location = None

    if log_scale:
        for axis in [fig.xaxis, fig.yaxis]:
            axis[0].ticker.base = log_scale
            axis[0].formatter.ticker = axis[0].ticker

    fig.grid.visible = grid
    fig.grid.name = 'grid'
    
    lasso = bokeh.models.LassoSelectTool(select_every_mousemove=False)
    fig.add_tools(lasso)
    
    if initial_xy_names is None:
        x_name, y_name = numerical_cols[:2]
    else:
        x_name, y_name = initial_xy_names
    
    fig.xaxis.name = 'x_axis'
    fig.yaxis.name = 'y_axis'
    fig.xaxis.axis_label = x_name
    fig.yaxis.axis_label = y_name
    for axis in (fig.xaxis, fig.yaxis):
        axis.axis_label_text_font_size = '{0}pt'.format(axis_label_size)
        axis.axis_label_text_font_style = 'normal'

    scatter_source.data['x'] = scatter_source.data[x_name]
    scatter_source.data['y'] = scatter_source.data[y_name]
    
    scatter_source.data['index'] = list(df.index)

    scatter_source.data['_no_color'] = ['rgba(0, 0, 0, 1.0)' for _ in scatter_source.data['x']]
    
    if color_by is None:
        color_by = '_no_color'
        show_color_by_menu = False
    else:
        show_color_by_menu = True

    if isinstance(color_by, basestring):
        color_options = ['', color_by]
    else:
        show_color_by_menu = True
        color_options = [''] + list(color_by)
    
    scatter_source.data['_color'] = scatter_source.data[color_options[1]]
    
    if label_by is None:
        label_by = df.index.name

    if isinstance(label_by, basestring):
        show_label_by_menu = False
        label_options = [label_by]
    else:
        show_label_by_menu = True
        label_options = list(label_by)
    
    scatter_source.data['_label'] = scatter_source.data[label_options[0]]


    scatter_source.selected = build_selected(initial_indices)

    scatter = fig.scatter('x',
                          'y',
                          source=scatter_source,
                          size=marker_size,
                          fill_color='_color',
                          fill_alpha=0.5,
                          line_color=None,
                          nonselection_color='_color',
                          nonselection_alpha=0.1,
                          selection_color='color',
                          selection_alpha=0.9,
                          name='scatter',
                         )
    
    if log_scale:
        nonzero = df[df > 0]

        overall_max = nonzero.max(numeric_only=True).max()
        overall_min = nonzero.min(numeric_only=True).min()
        
        initial = (overall_min * 0.1, overall_max * 10)
        bounds = (overall_min * 0.001, overall_max * 1000)
    
        def log(x):
            return np.log(x) / np.log(log_scale)

        bins = {}
        for name in numerical_cols:
            if identical_bins:
                left = overall_min
                right = overall_max
            else:
                left = nonzero[name].min()
                right = nonzero[name].max()

            bins[name] = list(np.logspace(log(left), log(right), num_bins))

    else:
        overall_max = df.max(numeric_only=True).max()
        overall_min = df.min(numeric_only=True).min()
        
        extent = overall_max - overall_min
        overhang = extent * 0.05
        max_overhang = extent * 0.5

        initial = (overall_min - overhang, overall_max + overhang)
        bounds = (overall_min - max_overhang, overall_max + max_overhang)
        
        bins = {}
        for name in numerical_cols:
            if identical_bins:
                left = overall_min
                right = overall_max
            else:
                left = df[name].min()
                right = df[name].max()

            bins[name] = list(np.linspace(left, right, num_bins))

    if data_lims is not None:
        initial = data_lims
        bounds = data_lims

    diagonals_visible = not grid # i.e. normally True

    fig.line(x=bounds, y=bounds,
             color='black',
             nonselection_color='black',
             alpha=0.4,
             nonselection_alpha=0.4,
             name='diagonal',
             visible=diagonals_visible,
            )

    if log_scale:
        upper_ys = np.array(bounds) * 10
        lower_ys = np.array(bounds) * 0.1
    else:
        upper_ys = np.array(bounds) + 1
        lower_ys = np.array(bounds) - 1

    line_kwargs = dict(
        color='black',
        nonselection_color='black',
        alpha=0.4,
        nonselection_alpha=0.4,
        line_dash=[5, 5],
        name='diagonal',
        visible=diagonals_visible,
    ) 
    fig.line(x=bounds, y=upper_ys, **line_kwargs)
    fig.line(x=bounds, y=lower_ys, **line_kwargs)
    
    if volcano:
        fig.y_range = bokeh.models.Range1d(-0.1, 8)
        fig.x_range = bokeh.models.Range1d(-1, 1)
    else:
        fig.y_range = bokeh.models.Range1d(*initial)
        fig.x_range = bokeh.models.Range1d(*initial)

    fig.x_range.name = 'x_range'
    fig.y_range.name = 'y_range'
    
    lower_bound, upper_bound = bounds
    range_kwargs = dict(lower_bound=lower_bound, upper_bound=upper_bound)
    
    fig.x_range.callback = build_callback('scatter_range', format_kwargs=range_kwargs)
    fig.y_range.callback = build_callback('scatter_range', format_kwargs=range_kwargs)
    
    fig.outline_line_color = 'black'

    scatter.selection_glyph.fill_color = 'orange'
    scatter.selection_glyph.line_color = None
    scatter.nonselection_glyph.line_color = None

    # Make marginal histograms.

    histogram_source = bokeh.models.ColumnDataSource(name='histogram_source')
    histogram_source.data = {
        'zero': [0]*(num_bins - 1),
    }

    for name in numerical_cols:
        histogram_source.data.update({
            '{0}_bins_left'.format(name): bins[name][:-1],
            '{0}_bins_right'.format(name): bins[name][1:],
        })

    max_count = 0
    for name in numerical_cols:
        counts, _ = np.histogram(df[name].dropna(), bins=bins[name])
        max_count = max(max(counts), max_count)
        histogram_source.data['{0}_all'.format(name)] = list(counts)

    if log_scale:
        axis_type = 'log'
    else:
        axis_type = 'linear'

    hist_figs = {
        'top': bokeh.plotting.figure(width=size, height=100,
                                     x_range=fig.x_range,
                                     x_axis_type=axis_type,
                                    ),
        'right': bokeh.plotting.figure(width=100, height=size,
                                       y_range=fig.y_range,
                                       y_axis_type=axis_type,
                                      ),
    }

    for axis, name in [('x', x_name), ('y', y_name)]:
        for data_type in ['all', 'bins_left', 'bins_right']:
            left_key = '{0}_{1}'.format(axis, data_type)
            right_key = '{0}_{1}'.format(name, data_type)
            histogram_source.data[left_key] = histogram_source.data[right_key]

        initial_vals = df[name].iloc[initial_indices]
        initial_counts, _ = np.histogram(initial_vals.dropna(), bins[name])
    
        histogram_source.data['{0}_selected'.format(axis)] = initial_counts

    initial_hist_alpha = 0.1 if len(initial_indices) > 0 else 0.2
    quads = {}
    quads['top_all'] = hist_figs['top'].quad(left='x_bins_left',
                                             right='x_bins_right',
                                             bottom='zero',
                                             top='x_all',
                                             source=histogram_source,
                                             color='black',
                                             alpha=initial_hist_alpha,
                                             line_color=None,
                                             name='hist_x_all',
                                            )

    quads['top_selected'] = hist_figs['top'].quad(left='x_bins_left',
                                                  right='x_bins_right',
                                                  bottom='zero',
                                                  top='x_selected',
                                                  source=histogram_source,
                                                  color='orange',
                                                  alpha=0.8,
                                                  line_color=None,
                                                 )

    quads['right_all'] = hist_figs['right'].quad(top='y_bins_left',
                                                 bottom='y_bins_right',
                                                 left='zero',
                                                 right='y_all',
                                                 source=histogram_source,
                                                 color='black',
                                                 alpha=initial_hist_alpha,
                                                 line_color=None,
                                                 name='hist_y_all',
                                                )

    quads['right_selected'] = hist_figs['right'].quad(top='y_bins_left',
                                                      bottom='y_bins_right',
                                                      left='zero',
                                                      right='y_selected',
                                                      source=histogram_source,
                                                      color='orange',
                                                      alpha=0.8,
                                                      line_color=None,
                                                     )

    # Poorly-understood bokeh behavior causes selection changes to sometimes
    # be broadcast to histogram_source. To prevent this from having any visual
    # effect, remove any difference in selection/nonselection glyphs.
    for quad in quads.values():
        quad.selection_glyph = quad.glyph
        quad.nonselection_glyph = quad.glyph

    hist_range_kwargs = dict(start=0, end=max_count, bounds='auto')
    hist_figs['top'].y_range = bokeh.models.Range1d(**hist_range_kwargs)
    hist_figs['right'].x_range = bokeh.models.Range1d(**hist_range_kwargs)

    for hist_fig in hist_figs.values():
        hist_fig.outline_line_color = None
        hist_fig.axis.visible = False
        hist_fig.grid.visible = False
        hist_fig.min_border = 0
        hist_fig.toolbar_location = None

    # Configure tooltips that pop up when hovering over a point.
    
    hover = bokeh.models.HoverTool()
    hover.tooltips = [
        (df.index.name, '@{0}'.format(df.index.name)),
    ]
    for key in hover_keys:
        hover.tooltips.append((key, '@{0}'.format(key)))
    fig.add_tools(hover)

    # Set up the table.

    table_col_names = [df.index.name] + table_keys
    columns = []
    for col_name in table_col_names:
        lengths = [len(str(v)) for v in scatter_source.data[col_name]]
        mean_length = np.mean(lengths)

        if col_name in numerical_cols:
            formatter = bokeh.models.widgets.NumberFormatter(format='0.00')
            width = 50
        else:
            formatter = None
            width = min(500, int(12 * mean_length))

        column = bokeh.models.widgets.TableColumn(field=col_name,
                                                  title=col_name,
                                                  formatter=formatter,
                                                  width=width,
                                                 )
        columns.append(column)

    filtered_data = {k: [scatter_source.data[k][i] for i in initial_indices]
                     for k in scatter_source.data
                    }
    
    filtered_source = bokeh.models.ColumnDataSource(data=filtered_data, name='labels_source')
    
    table = bokeh.models.widgets.DataTable(source=filtered_source,
                                           columns=columns,
                                           width=2 * size if heatmap else size,
                                           height=600,
                                           sortable=False,
                                           name='table',
                                           row_headers=False,
                                          )
    
    # Callback to filter the table when selection changes.
    scatter_source.callback = build_callback('scatter_selection')
    
    labels = bokeh.models.LabelSet(x='x',
                                   y='y',
                                   text='_label',
                                   level='glyph',
                                   x_offset=0,
                                   y_offset=2,
                                   source=filtered_source,
                                   text_font_size='8pt',
                                   name='labels',
                                  )
    fig.add_layout(labels)
    
    # Set up menus or heatmap to select columns from df to put on x- and y-axis.

    if heatmap:
        norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)

        def r_to_color(r):
            color = matplotlib.colors.rgb2hex(matplotlib.cm.RdBu_r(norm(r)))
            return color

        data = {
            'x': [],
            'x_name': [],
            'y': [],
            'y_name': [],
            'r': [],
            'color': [],
        }

        correlations = df.corr()
        for y, row in enumerate(numerical_cols):
            for x, col in enumerate(numerical_cols):
                r = correlations[row][col]
                data['r'].append(r)
                data['x'].append(x)
                data['x_name'].append(col)
                data['y'].append(y)
                data['y_name'].append(row)
                data['color'].append(r_to_color(r))
                
        heatmap_source = bokeh.models.ColumnDataSource(data)
        num_exps = len(numerical_cols)
        heatmap_size = size - 100
        heatmap_fig = bokeh.plotting.figure(tools='tap',
                                            x_range=(-0.5, num_exps - 0.5),
                                            y_range=(num_exps - 0.5, -0.5),
                                            width=heatmap_size, height=heatmap_size,
                                            toolbar_location=None,
                                           )

        heatmap_fig.grid.visible = False
        rects = heatmap_fig.rect(x='x', y='y',
                                 line_color=None,
                                 hover_line_color='black',
                                 hover_fill_color='color',
                                 selection_fill_color='color',
                                 nonselection_fill_color='color',
                                 nonselection_fill_alpha=1,
                                 nonselection_line_color=None,
                                 selection_line_color='black',
                                 line_width=5,
                                 fill_color='color',
                                 source=heatmap_source,
                                 width=1, height=1,
                                )

        hover = bokeh.models.HoverTool()
        hover.tooltips = [
            ('X', '@x_name'),
            ('Y', '@y_name'),
            ('r', '@r'),
        ]
        heatmap_fig.add_tools(hover)

        first_row = [heatmap_fig]
        heatmap_source.callback = build_callback('scatter_heatmap')

        code = '''
        dict = {dict}
        return dict[tick].slice(0, 15)
        '''.format(dict=dict(enumerate(numerical_cols)))
        
        for ax in [heatmap_fig.xaxis, heatmap_fig.yaxis]:
            ax.ticker = bokeh.models.FixedTicker(ticks=range(num_exps))
            ax.formatter = bokeh.models.FuncTickFormatter(code=code)
            ax.major_label_text_font_size = '8pt'

        heatmap_fig.xaxis.major_label_orientation = np.pi / 4

        name_pairs = zip(heatmap_source.data['x_name'], heatmap_source.data['y_name'])
        initial_index = name_pairs.index((x_name, y_name))
        heatmap_source.selected = build_selected([initial_index])

        heatmap_fig.min_border = min_border
        
    else:
        x_menu = bokeh.models.widgets.MultiSelect(title='X',
                                                  options=numerical_cols,
                                                  value=[x_name],
                                                  size=min(6, len(numerical_cols)),
                                                  name='x_menu',
                                               )
        y_menu = bokeh.models.widgets.MultiSelect(title='Y',
                                                  options=numerical_cols,
                                                  value=[y_name],
                                                  size=min(6, len(numerical_cols)),
                                                  name='y_menu',
                                               )

        menu_callback = build_callback('scatter_menu')
        x_menu.callback = menu_callback
        y_menu.callback = menu_callback
        
        first_row = [bokeh.layouts.widgetbox([x_menu, y_menu])],
    
    # Button to toggle labels.
    label_button = bokeh.models.widgets.Toggle(label='label selected points',
                                               width=50,
                                               active=True,
                                               name='label_button',
                                              )
    label_button.callback = bokeh.models.CustomJS(args={'labels': labels},
                                                  code='labels.text_alpha = 1 - labels.text_alpha;',
                                                 )
    
    # Button to zoom to current data limits.
    zoom_to_data_button = bokeh.models.widgets.Button(label='zoom to data limits',
                                                      width=50,
                                                      name='zoom_button',
                                                     )
    zoom_to_data_button.callback = build_callback('scatter_zoom_to_data',
                                                  format_kwargs=dict(log_scale='true' if log_scale else 'false'),
                                                 )
    
    
    # Menu to choose label source.
    label_menu = bokeh.models.widgets.Select(title='Label by:',
                                             options=label_options,
                                             value=label_options[0],
                                             name='label_menu',
                                            )
    label_menu.callback = build_callback('scatter_label')

    # Menu to choose color source.
    color_menu = bokeh.models.widgets.Select(title='Color by:',
                                             options=color_options,
                                             value=color_options[1],
                                             name='color_menu',
                                            )
    color_menu.callback = build_callback('scatter_color')

    # Radio group to choose whether to draw a vertical/horizontal grid or
    # diagonal guide lines. 
    grid_options = bokeh.models.widgets.RadioGroup(labels=['grid', 'diagonal'],
                                                   active=1 if not grid else 0,
                                                   name='grid_radio_buttons',
                                                  )
    grid_options.callback = build_callback('scatter_grid')

    text_input = bokeh.models.widgets.TextInput(title='Search:', name='search')
    text_input.callback = build_callback('scatter_search',
                                                format_kwargs=dict(column_names=str(object_cols)),
                                               )

    case_sensitive = bokeh.models.widgets.CheckboxGroup(labels=['Case sensitive'],
                                                        active=[],
                                                        name='case_sensitive',
                                                       )
    case_sensitive.callback = build_callback('case_sensitive')

    # Menu to select a subset of points from a columns of bools.
    subset_options = [''] + bool_cols
    subset_menu = bokeh.models.widgets.Select(title='Select subset:',
                                              options=subset_options,
                                              value='',
                                              name='subset_menu',
                                             )
    subset_menu.callback = build_callback('scatter_subset_menu',
                                          format_kwargs=dict(subset_indices=str(subset_indices)),
                                         )

    # Button to dump table to file.
    save_button = bokeh.models.widgets.Button(label='Save table to file',
                                              width=50,
                                              name='save_button',
                                             )
    save_button.callback = build_callback('scatter_save_button',
                                          format_kwargs=dict(column_names=str(table_col_names)),
                                         )

    if alpha_widget_type == 'slider':
        alpha_widget = bokeh.models.Slider(start=0.,
                                           end=1.,
                                           value=0.5,
                                           step=.05,
                                           title='alpha',
                                           name='alpha',
                                          )
    elif alpha_widget_type == 'text':
        alpha_widget = bokeh.models.TextInput(title='alpha', name='alpha', value='0.5')
    else:
        raise ValueError('{0} not a valid alpha_widget_type value'.format(alpha_widget_type))

    alpha_widget.callback = build_callback('scatter_alpha')
    
    size_slider = bokeh.models.Slider(start=1,
                                      end=20.,
                                      value=marker_size,
                                      step=1,
                                      title='marker size',
                                      name='marker_size',
                                      )
    size_slider.callback = build_callback('scatter_size')

    fig.min_border = 1

    widgets = [
        label_button,
        zoom_to_data_button,
        label_menu,
        color_menu,
        grid_options,
        alpha_widget,
        size_slider,
        text_input,
        case_sensitive,
        subset_menu,
        save_button,
    ]

    hide_widgets = set(hide_widgets)

    if 'table' in hide_widgets:
        hide_widgets.add('save_button')

    if 'search' in hide_widgets:
        hide_widgets.add('case_sensitive')

    if not show_color_by_menu:
        hide_widgets.add('color_menu')

    if not show_label_by_menu:
        hide_widgets.add('label_menu')

    if len(subset_options) == 1:
        hide_widgets.add('subset_menu')

    widgets = [w for w in widgets if w.name not in hide_widgets]

    if not heatmap:
        widgets = [x_menu, y_menu] + widgets

    widget_box = bokeh.layouts.widgetbox(children=widgets)

    toolbar = bokeh.models.ToolbarBox(tools=fig.toolbar.tools)
    toolbar.logo = None

    columns = [
        bokeh.layouts.column(children=[hist_figs['top'], fig]),
        bokeh.layouts.column(children=[bokeh.layouts.Spacer(height=100), hist_figs['right']]),
        bokeh.layouts.column(children=[bokeh.layouts.Spacer(height=100), toolbar]),
        bokeh.layouts.column(children=[bokeh.layouts.Spacer(height=min_border),
                                       widget_box,
                                      ]),
    ]

    if heatmap:
        heatmap_column = bokeh.layouts.column(children=[bokeh.layouts.Spacer(height=100), heatmap_fig])
        columns = columns[:-1] + [heatmap_column] + columns[-1:]

    rows = [
        bokeh.layouts.row(children=columns),
    ]
    if 'table' not in hide_widgets:
        rows.append(table)

    full = bokeh.layouts.column(children=rows)

    bokeh.io.show(full)

def hex_to_CSS(hex_string, alpha=1.):
    ''' Converts an RGB hex value and optional alpha value to a CSS-format RGBA string. '''
    rgb = matplotlib.colors.colorConverter.to_rgb(hex_string)
    rgb = [int(v * 255) for v in rgb]
    CSS = 'rgba({1}, {2}, {3}, {0})'.format(alpha, *rgb)
    return CSS

# from http://chris-said.io/
toggle = '''
<script>
  function code_toggle() {
    if (code_shown){
      $('div.input').hide('100');
      $('#toggleButton').val('Show Code')
    } else {
      $('div.input').show('100');
      $('#toggleButton').val('Hide Code')
    }
    code_shown = !code_shown
  }

  $( document ).ready(function(){
    code_shown=false;
    $('div.input').hide()
  });
</script>
<form action="javascript:code_toggle()"><input type="submit" id="toggleButton" value="Show Code"></form>
'''

toggle_cell = IPython.display.HTML(toggle)
