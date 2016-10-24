import numpy as np
import bokeh.io
import bokeh.plotting
import pandas as pd
import PIL.ImageColor

bokeh.io.output_notebook()

def scatter(df, hover_keys=None, table_keys=None, size=900):
    ''' Makes an interactive scatter plot using bokeh.

    Args:
            df: A pandas DataFrame with columns containing numerical data to plot.
                If 'color' is a column, it will be used to color the points. 
                Index values will be used as labels for points.
            hover_keys: Names of columns in df to display in the tooltip that appears
                when you hover over a point.
            table_keys: Names of columns in df to display in the table below the plot
                that is populated with the selected points from the figure.
    '''

    if hover_keys is None:
        hover_keys = []

    if table_keys is None:
        table_keys = []


    # Set up the actual scatter plot.
    
    tools= [
        'reset',
        'pan',
        'box_zoom',
        'box_select',
        'poly_select',
        'tap',
        'wheel_zoom',
        'save',
    ]
    
    fig = bokeh.plotting.figure(plot_width=size,
                                plot_height=size,
                                tools=','.join(tools),
                                lod_threshold=5000,
                               )
    
    lasso = bokeh.models.LassoSelectTool(select_every_mousemove=False)
    fig.add_tools(lasso)

    numerical_cols = [n for n in df.columns if df[n].dtype in [float, int]]
    x_name, y_name = numerical_cols[:2]
    
    fig.xaxis.axis_label = x_name
    fig.yaxis.axis_label = y_name
    for axis in (fig.xaxis, fig.yaxis):
        axis.axis_label_text_font_size = '20pt'
        axis.axis_label_text_font_style = 'normal'

    scatter_source = bokeh.models.ColumnDataSource(data=df)

    scatter_source.data['x'] = scatter_source.data[x_name]
    scatter_source.data['y'] = scatter_source.data[y_name]

    if 'color' not in df:
        scatter_source.data['color'] = ['rgba(0, 0, 0, 0.5)' for _ in scatter_source.data['x']]

    if df.index.name is None:
        df.index.name = 'index'
    
    scatter = fig.scatter('x',
                          'y',
                          source=scatter_source,
                          size=6,
                          fill_color='color',
                          line_color=None,
                         )
    
    overall_max = df.max(numeric_only=True).max()
    overall_min = df.min(numeric_only=True).min()
    
    extent = overall_max - overall_min
    overhang = extent * 0.05
    max_overhang = extent * 0.5
    
    initial = (overall_min - overhang, overall_max + overhang)
    bounds = (overall_min - max_overhang, overall_max + max_overhang)
    
    fig.y_range = bokeh.models.Range1d(*initial, bounds=bounds)
    fig.x_range = bokeh.models.Range1d(*initial, bounds=bounds)

    if 'color' not in df:
        scatter.selection_glyph = bokeh.models.Circle(fill_color="orange",
                                                      line_color=None,
                                                     )

    scatter.nonselection_glyph = bokeh.models.Circle(fill_color="black",
                                                     line_color=None,
                                                     fill_alpha=0.1,
                                                    )

    # Configure tooltips that pop up when hovering over a point.
    
    hover = bokeh.models.HoverTool()
    hover.tooltips = [
        (df.index.name, '@{0}'.format(df.index.name)),
        ('x', '@x'),
        ('y', '@y'),
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
    filtered_source = bokeh.models.ColumnDataSource(data=filtered_data)
    
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
                                           width=size,
                                           height=1000,
                                           sortable=False,
                                          )
    
    # Set up menus to select columns from df to put on x- and y-axis.

    x_menu = bokeh.models.widgets.Select(title='X', options=numerical_cols, value=x_name)
    y_menu = bokeh.models.widgets.Select(title='Y', options=numerical_cols, value=y_name)

    menu_callback_code = """
    var scatter_data = scatter_source.data;
    var label_data = label_source.data;
    
    var x_name = x_menu.value;
    var y_name = y_menu.value;
    
    scatter_data.x = scatter_data[x_name];
    scatter_data.y = scatter_data[y_name];
    
    label_data.x = label_data[x_name];
    label_data.y = label_data[y_name];
    
    xaxis.axis_label = x_name;
    yaxis.axis_label = y_name;
    
    scatter_source.trigger('change');
    label_source.trigger('change');
    """

    menu_args = dict(scatter_source=scatter.data_source,
                     label_source=labels.source,
                     x_menu=x_menu,
                     y_menu=y_menu,
                     xaxis=fig.xaxis[0],
                     yaxis=fig.yaxis[0],
                    )
    menu_callback = bokeh.models.CustomJS(args=menu_args, code=menu_callback_code)
    x_menu.callback = menu_callback
    y_menu.callback = menu_callback
    
    # Set up callback to filter the table when selection changes.

    selection_callback_code = """
    full_data = source.data
    filtered_data = table.source.data
    indices = cb_obj.selected['1d'].indices

    for key, values of full_data
        filtered_data[key] = (values[i] for i in indices)

    table.trigger('change')
    labels.trigger('change')
    """

    selection_args = dict(source=scatter_source, table=table, labels=labels)
    scatter_source.callback = bokeh.models.CustomJS.from_coffeescript(args=selection_args, code=selection_callback_code)
    
    # Button to toggle labels.
    
    button = bokeh.models.widgets.Toggle(label='label selected points',
                                         width=50,
                                         active=True,
                                        )
    button.callback = bokeh.models.CustomJS(args={'labels': labels},
                                            code='labels.text_alpha = 1 - labels.text_alpha;',
                                           )

    grid = [
        [x_menu, y_menu, button],
        [fig],
        [table],
    ]
    layout = bokeh.layouts.layout(grid)
    bokeh.io.show(layout)

def hex_to_CSS(hex_string, alpha=1.):
    ''' Converts an RGB hex value and option alpha value to a CSS-format RGBA string. '''
    rgb = PIL.ImageColor.getrgb(hex_string)
    CSS = 'rgba({1}, {2}, {3}, {0})'.format(alpha, *rgb)
    return CSS

def example():
    df = pd.read_csv('/home/jah/projects/sequencing/code/Sequencing/Visualize/example_df.txt', index_col='alias')
    scatter(df, hover_keys=['short_description'], table_keys=['description'])
