import numpy as np
import bokeh.io
import bokeh.plotting
import pandas as pd
import scipy.stats
import matplotlib.colors
import matplotlib.cm
import os.path
import glob
from collections import defaultdict

bokeh.io.output_notebook()

# For easier editing, coffeescript callbacks are kept in separate files
# in the same directory as this one. Load their contents into a dictionary.

coffee_fns = glob.glob(os.path.join(os.path.dirname(__file__), '*.coffee'))
callbacks = {}
for fn in coffee_fns:
    head, tail = os.path.split(fn)
    root, ext = os.path.splitext(tail)
    with open(fn) as fh:
        callbacks[root] = fh.read()

def external_coffeescript(key, args=None, format_args=None):
    if args is None:
        args = {}
    if format_args is None:
        format_args = {}

    code = callbacks[key].format(**format_args)
    callback = bokeh.models.CustomJS.from_coffeescript(code=code, args=args)
    return callback

colors_list =  (
    bokeh.palettes.Dark2[8] +
    bokeh.palettes.Set1[9] +
    bokeh.palettes.Set2[8] + 
    bokeh.palettes.Paired[12] +
    bokeh.palettes.Accent[8]
)

def scatter(df,
            hover_keys=None, table_keys=None,
            size=900,
            axis_label_size=20,
            log_scale=False,
            volcano=False,
            heatmap=False,
            grid=False,
            marker_size=6,
           ):
    ''' Makes an interactive scatter plot using bokeh.

    Args:
            df: A pandas DataFrame with columns containing numerical data to plot.
                If 'color' is a column, it will be used to color the points. 
                Index values will be used as labels for points.
                Any text columns will be searchable through the 'Search:' field.
                Any boolean columns will be used to define subset of points for
                selection from a dropdown menu.
            hover_keys: Names of columns in df to display in the tooltip that appears
                when you hover over a point.
            table_keys: Names of columns in df to display in the table below the plot
                that is populated with the selected points from the figure.
            size: Size of the plot in pixels.
            marker_size: Size of the scatter circles.
            heatmap: If True, displays a heatmap of correlations between
                numerical columns in df that can be clicked to select columns
                to scatter.
            grid: If True, defaults to grid instead of diagonal landmarks.
            volcano: If True, make some tweaks suitable for volcano plots.
            log_scale: If not False, plot on a log scale with base 10 (or, if a
                set to a number, with base log_scale.)
    '''

    if hover_keys is None:
        hover_keys = []

    if table_keys is None:
        table_keys = []

    if volcano:
        grid = True

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

    if log_scale:
        if log_scale is True:
            log_scale = 10
        fig_kwargs['y_axis_type'] = 'log'
        fig_kwargs['x_axis_type'] = 'log'
    
    fig = bokeh.plotting.figure(**fig_kwargs)
    fig.toolbar.logo = None

    if log_scale:
        for axis in [fig.xaxis, fig.yaxis]:
            axis[0].ticker.base = log_scale
            axis[0].formatter.ticker = axis[0].ticker

    fig.grid.visible = grid
    fig.grid.name = 'grid'
    
    lasso = bokeh.models.LassoSelectTool(select_every_mousemove=False)
    fig.add_tools(lasso)

    numerical_cols = [n for n in df.columns if df[n].dtype in [float, int]]

    object_cols = [n for n in df.columns if df[n].dtype is np.dtype('O')]
    if df.index.dtype is np.dtype('O'):
        object_cols.append(df.index.name)

    bool_cols = [n for n in df.columns if df[n].dtype is np.dtype('bool')]

    x_name, y_name = numerical_cols[:2]
    
    fig.xaxis.axis_label = x_name
    fig.yaxis.axis_label = y_name
    for axis in (fig.xaxis, fig.yaxis):
        axis.axis_label_text_font_size = '{0}pt'.format(axis_label_size)
        axis.axis_label_text_font_style = 'normal'

    scatter_source = bokeh.models.ColumnDataSource(data=df, name='scatter_source')

    scatter_source.data['x'] = scatter_source.data[x_name]
    scatter_source.data['y'] = scatter_source.data[y_name]

    scatter_source.data['no_color'] = ['rgba(0, 0, 0, 0.5)' for _ in scatter_source.data['x']]
    if 'color' not in scatter_source.data:
        scatter_source.data['color'] = scatter_source.data['no_color']

    if df.index.name is None:
        df.index.name = 'index'
    
    scatter = fig.scatter('x',
                          'y',
                          source=scatter_source,
                          size=marker_size,
                          fill_color='color',
                          line_color=None,
                          nonselection_color='color',
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
    else:
        overall_max = df.max(numeric_only=True).max()
        overall_min = df.min(numeric_only=True).min()

        extent = overall_max - overall_min
        overhang = extent * 0.05
        max_overhang = extent * 0.5

        initial = (overall_min - overhang, overall_max + overhang)
        bounds = (overall_min - max_overhang, overall_max + max_overhang)

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
        fig.y_range = bokeh.models.Range1d(*initial, bounds=bounds)
        fig.x_range = bokeh.models.Range1d(*initial, bounds=bounds)

    fig.x_range.name = 'x_range'
    fig.y_range.name = 'y_range'

    scatter.selection_glyph.fill_color = "orange"
    scatter.selection_glyph.line_color = None
    scatter.nonselection_glyph.line_color = None

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
        if col_name == df.index.name:
            formatter = None
            width = 80
        elif col_name in numerical_cols:
            formatter = bokeh.models.widgets.NumberFormatter(format='0.00')
            width = 50
        else:
            formatter = None
            width = 500

        column = bokeh.models.widgets.TableColumn(field=col_name,
                                                  title=col_name,
                                                  formatter=formatter,
                                                  width=width,
                                                 )
        columns.append(column)

    filtered_data = {k: [] for k in list(df.columns) + [df.index.name, 'x', 'y']}
    filtered_source = bokeh.models.ColumnDataSource(data=filtered_data, name='labels_source')
    
    labels = bokeh.models.LabelSet(x='x',
                                   y='y',
                                   text=df.index.name,
                                   level='glyph',
                                   x_offset=0,
                                   y_offset=2,
                                   source=filtered_source,
                                   text_font_size='8pt',
                                  )
    fig.add_layout(labels)
    
    table = bokeh.models.widgets.DataTable(source=filtered_source,
                                           columns=columns,
                                           width=2 * size if heatmap else size,
                                           height=1000,
                                           sortable=False,
                                           name='table',
                                          )
    
    # Set up menus or heatmap to select columns from df to put on x- and y-axis.

    if heatmap:
        norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)
        def r_to_color(r):
            color = matplotlib.colors.rgb2hex(matplotlib.cm.RdBu_r(norm(r)))
            return color
        xs = []
        x_names = []
        ys = []
        y_names = []
        rs = []
        colors = []

        for y, row in enumerate(numerical_cols):
            for x, col in enumerate(numerical_cols):
                r, p = scipy.stats.pearsonr(df[row], df[col])
                rs.append(r)
                xs.append(x)
                x_names.append(col)
                ys.append(y)
                y_names.append(row)
                colors.append(r_to_color(r))
                
        data = {
            'x': xs,
            'x_name': x_names,
            'y': ys,
            'y_name': y_names,
            'r': rs,
            'color': colors,
        }

        heatmap_source = bokeh.models.ColumnDataSource(data)
        num_exps = len(numerical_cols)
        heatmap_size = int(size * 1)
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
        args = dict(xaxis=fig.xaxis[0],
                    yaxis=fig.yaxis[0],
                   )
        heatmap_source.callback = external_coffeescript('scatter_heatmap', args=args)

        code = '''
        dict = {dict}
        return dict[tick].slice(0, 15)
        '''.format(dict=dict(enumerate(numerical_cols)))
        
        for ax in [heatmap_fig.xaxis, heatmap_fig.yaxis]:
            ax.ticker = bokeh.models.FixedTicker(ticks=range(num_exps))
            ax.formatter = bokeh.models.FuncTickFormatter(code=code)
            ax.major_label_text_font_size = '8pt'

        heatmap_fig.xaxis.major_label_orientation = np.pi / 4

        pvd = bokeh.core.property_containers.PropertyValueDict
        pvl = bokeh.core.property_containers.PropertyValueList

        selected_from_scratch = pvd({
            '0d': pvd({
                'glyph': None,
                'indices': pvl(),
            }),
            '1d': pvd({
                'indices': pvl([num_exps]),
            }),
            '2d': pvd(),
        })

        heatmap_source.selected = selected_from_scratch
        heatmap_fig.min_border = 80
        
    else:
        x_menu = bokeh.models.widgets.Select(title='X',
                                             options=numerical_cols,
                                             value=x_name,
                                             name='x_menu',
                                            )
        y_menu = bokeh.models.widgets.Select(title='Y',
                                             options=numerical_cols,
                                             value=y_name,
                                             name='y_menu',
                                            )

        menu_args = dict(xaxis=fig.xaxis[0],
                         yaxis=fig.yaxis[0],
                        )
        menu_callback = external_coffeescript('scatter_menu', args=menu_args)
        x_menu.callback = menu_callback
        y_menu.callback = menu_callback
        
        first_row = [bokeh.layouts.widgetbox([x_menu, y_menu])],
    
    # Set up callback to filter the table when selection changes.
    scatter_source.callback = external_coffeescript('scatter_selection')
    
    # Button to toggle labels.
    label_button = bokeh.models.widgets.Toggle(label='label selected points',
                                               width=50,
                                               active=True,
                                              )
    label_button.callback = bokeh.models.CustomJS(args={'labels': labels},
                                                  code='labels.text_alpha = 1 - labels.text_alpha;',
                                                 )
    
    # Button to zoom to current data limits.
    zoom_to_data_button = bokeh.models.widgets.Button(label='zoom to data limits',
                                                      width=50,
                                                     )
    zoom_to_data_button.callback = external_coffeescript('scatter_zoom_to_data',
                                                         format_args=dict(log_scale='true' if log_scale else 'false'))

    grid_options = bokeh.models.widgets.RadioGroup(labels=['grid', 'diagonal'], active=1 if not grid else 0)
    grid_options.callback = external_coffeescript('scatter_grid')

    text_input = bokeh.models.widgets.TextInput(title='Search:')
    text_input.callback = external_coffeescript('scatter_search',
                                                format_args=dict(columns=object_cols),
                                               )
    
    # Menu to select a subset of points from a columns of bools.
    subset_options = [''] + bool_cols
    subset_menu = bokeh.models.widgets.Select(title='Select subset:',
                                             options=subset_options,
                                             value='',
                                            )
    subset_menu.callback = external_coffeescript('scatter_subset_menu')

    # Button to dump table to file.
    save_button = bokeh.models.widgets.Button(label='Save table to file',
                                              width=50,
                                             )
    save_button.callback = external_coffeescript('scatter_save_button',
                                                 format_args=dict(columns=table_col_names),
                                                )

    fig.min_border = 80

    widgets = [
        label_button,
        zoom_to_data_button,
        grid_options,
        text_input,
        subset_menu,
        save_button,
    ]

    if len(subset_options) == 1:
        widgets = widgets[:-2] + widgets[-1:]

    if not heatmap:
        widgets = [x_menu, y_menu] + widgets

    widget_box = bokeh.layouts.widgetbox(children=widgets)

    columns = [
        bokeh.layouts.column(children=[fig]),
        bokeh.layouts.column(children=[bokeh.layouts.Spacer(height=100),
                                       widget_box,
                                      ]),
    ]

    if heatmap:
        columns = columns[:1] + [bokeh.layouts.column(children=[heatmap_fig])] + columns[1:]

    full = bokeh.layouts.column(children=[
        bokeh.layouts.row(children=columns),
        table,
    ])

    bokeh.io.show(full)

def hex_to_CSS(hex_string, alpha=1.):
    ''' Converts an RGB hex value and option alpha value to a CSS-format RGBA string. '''
    rgb = matplotlib.colors.colorConverter.to_rgb(hex_string)
    CSS = 'rgba({1}, {2}, {3}, {0})'.format(alpha, *rgb)
    return CSS

def example():
    fn = os.path.join(os.path.dirname(__file__), 'example_df.txt') 
    df = pd.read_csv(fn, index_col='alias')
    scatter(df, heatmap=True, hover_keys=['short_description'], table_keys=['description'], grid=True)
