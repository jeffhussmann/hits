import numbers
import base64
import re
from collections import defaultdict
from pathlib import Path

from matplotlib.pyplot import hist

import numpy as np
import scipy.cluster.hierarchy
import bokeh
import bokeh.io
import bokeh.plotting
import pandas as pd
import matplotlib.colors
import matplotlib.cm
import IPython.display

bokeh.io.output_notebook()

def bool_to_js(b):
    return 'true' if b else 'false'

from ..callback import build_js_callback

def scatter(df=None,
            numerical_cols=None,
            hover_keys=None,
            table_keys=None,
            color_by=None,
            label_by=None,
            size=800,
            label_size=8,
            axis_label_size=20,
            log_scale=False,
            volcano=False,
            heatmap=False,
            cluster=False,
            grid='none',
            marker_size=6,
            initial_selection=None,
            initial_xy_names=None,
            initial_alpha=0.5,
            nonselection_alpha=0.1,
            initial_data_lims=None,
            data_bounds=None,
            zoom_to_initial_data=False,
            alpha_widget_type='slider',
            hide_widgets=None,
            identical_bins=False,
            num_bins=100,
            show_axes_lines=False,
            return_layout=False,
            title=None,
            two_level_index=False,
            menu_width=200,
            save_as=None,
            level_order=None,
           ):
    ''' Makes an interactive scatter plot using bokeh. Call without any
    arguments for an example using data from Jan et al. Science 2014.

    Args:
        df: A pandas DataFrame with columns containing numerical data to plot. 
            Index values will be used as labels for points.
            Any text columns will be searchable through the 'Search:' field.
            Any boolean columns will be used to define subsets of points for
            selection from a dropdown menu.
            If df is None, loads example data from Jan et al. Science 2014.

        numerical_cols: If given, a list of columns to use as plotting choices.
            (If not given, all columns containing numerical data will be used.)

        hover_keys: Names of columns in df to display in the tooltip that
            appears when you hover over a point.

        table_keys: Names of columns in df to display in the table below the
            plot that is populated with the selected points from the figure.

        color_by: The name of a column in df to use as colors of points, or a
            list of such names to choose from a menu. (These columns will be
            excluded from text searches.)
        
        label_by: The name of a column in df to use as labels of points, or a
            list of such names to choose from a menu. If None, df.index is used.

        size: Size of the plot in pixels.

        marker_size: Size of the scatter circles. Can be a scalar value, a
            column name, or a list of column names to choose from via dropdown
            menu.

        label_size: Size in pts of font for labels of selected point.

        heatmap: If True, displays a heatmap of correlations between numerical
            columns in df that can be clicked to select columns to scatter. If
            any negative correlations exist, uses a red/blue colormap, otherwise
            uses viridis.

        cluster: If True, forces heatmap=True and performs hierarchical
            clustering on (correlations of) numerical columns and draws a
            dendrogram.

        grid: Draw a 'grid', 'diagonal' lines, 'diagonal+axes', or 'none' as guide lines.

        volcano: If True, make some tweaks suitable for volcano plots.

        log_scale: If not False, plot on a log scale with base 10 (or, if a set
            to a number, with base log_scale.)

        axis_label_size: Size of the font used for axis labels.

        initial_selection: Names of index value to initially highlight.

        initial_xy_names: Tuple (x_name, y_name) of datasets to initially
            display on x- and y-axes.

        alpha_widget_type: Type of widget ('slider' or 'text') to control
            scatter circle transparency.

        hide_widgets: List of widgets to not display. Possible options are
            ['table', 'alpha', 'marker_size', 'search', 'subset_menu',
             'grid_radio_buttons'].

        identical_bins: If True, use the same set of bins for histograms of all
            data sets. If False, use bins specific to the range of each data
            set.

        num_bins: Number of bins to use for marginal histograms.

        zoom_to_initial_data: If True, zoom to data limits of the initially
            selected columns rather than global data limits.

        return_layout: If True, return the final layout object to allowing
            embedding.
    '''

    if save_as is not None:
        bokeh.io.output_file(str(save_as))

    if hover_keys is None:
        hover_keys = []

    if table_keys is None:
        table_keys = []

    if level_order is None:
        level_order = {}

    if hide_widgets is None:
        hide_widgets = set()
    else:
        hide_widgets = set(hide_widgets)

    if volcano:
        grid = 'grid'

    if cluster:
        heatmap = True

    if df is None:
        # Load example data.
        fn = Path(__file__).parent / 'example_df.txt'
        df = pd.read_csv(fn, index_col='alias')

        # Override some arguments.
        log_scale = True
        hover_keys = ['systematic_name', 'short_description']
        table_keys = ['systematic_name', 'description']
        grid = 'diagonal'
        heatmap = False
        identical_bins = True

    # Copy before changing.
    original_index_name = df.index.name
    original_df = df
    df = df.copy()

    # Collapse multiindex if present
    df.columns = [' '.join(map(str, n)) if isinstance(n, tuple) else n for n in df.columns]
    df.index = [' '.join(map(str, n)) if isinstance(n, tuple) else n for n in df.index.values]
    df.index.name = original_index_name

    # Infer column types.
    scatter_data = df.to_dict(orient='list')

    if 'index' in scatter_data:
        scatter_data['_index'] = scatter_data['index']

    if df.index.name is None:
        df.index.name = 'index'
    
    scatter_data[df.index.name] = list(df.index)

    if initial_selection is None:
        initial_selection = []

    initial_indices = [i for i, n in enumerate(df.index) if n in initial_selection]

    auto_numerical_cols = list(df.select_dtypes('number').columns)
    if numerical_cols is not None:
        for col in numerical_cols:
            if col not in auto_numerical_cols:
                raise ValueError(col + ' not a numerical column')
    else:
        numerical_cols = auto_numerical_cols

    # bokeh can handle NaNs in numpy arrays but not in lists. 
    for numerical_col in numerical_cols:
        scatter_data[numerical_col] = np.array(scatter_data[numerical_col])

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

    plot_width = size + (35 if two_level_index else 0) # empirical adjustment for extra blank text row
    plot_height = size
    
    fig_kwargs = dict(
        plot_width=plot_width,
        plot_height=plot_height,
        tools=tools,
        lod_threshold=40000,
        name='scatter_fig',
        toolbar_location=None,
    )

    min_border = 80

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

    fig.grid.visible = (grid == 'grid')
    fig.grid.name = 'grid'
    
    lasso = bokeh.models.LassoSelectTool(select_every_mousemove=False)
    fig.add_tools(lasso)
    
    if initial_xy_names is None:
        if two_level_index:
            first_level = original_df.columns.levels[0]
            second_level = [s for s in original_df.columns.levels[1] if s != '']

            if first_level[1] != 'gene':
                initial_xy_names = (
                    (first_level[0], second_level[0]),
                    (first_level[1], second_level[0]),
                )
            else:
                initial_xy_names = (
                    (first_level[0], second_level[0]),
                    (first_level[0], second_level[1]),
                )
        else:
            initial_xy_names = numerical_cols[:2]

    if two_level_index:
        x_name, y_name = [' '.join(t) for t in initial_xy_names]
    else:
        x_name, y_name = initial_xy_names
    
    fig.xaxis.name = 'x_axis'
    fig.yaxis.name = 'y_axis'

    if two_level_index:
        x_title = bokeh.models.Title(text=initial_xy_names[0][0], text_font_size='16pt', align='center')
        x_subtitle = bokeh.models.Title(text=initial_xy_names[0][1], text_font_size='14pt', align='center')
        fig.add_layout(x_title, 'below')
        fig.add_layout(x_subtitle, 'below')

        y_title = bokeh.models.Title(text=initial_xy_names[1][0], text_font_size='16pt', align='center')
        y_subtitle = bokeh.models.Title(text=initial_xy_names[1][1], text_font_size='14pt', align='center')

        fig.add_layout(bokeh.models.Title(text=' ', text_font_size='6pt'), 'left')
        fig.add_layout(y_subtitle, 'left')
        fig.add_layout(y_title, 'left')
    else:
        fig.xaxis.axis_label = x_name
        fig.yaxis.axis_label = y_name
        for axis in (fig.xaxis, fig.yaxis):
            axis.axis_label_text_font_size = f'{axis_label_size}pt'
            axis.axis_label_text_font_style = 'normal'

    scatter_data['x'] = scatter_data[x_name]
    scatter_data['y'] = scatter_data[y_name]
    
    scatter_data['index'] = list(df.index)

    scatter_data['_black'] = ['rgba(0, 0, 0, 1.0)' for _ in scatter_data['x']]
    scatter_data['_orange'] = ['orange' for _ in scatter_data['x']]
    
    if color_by is None:
        color_by = ''
        show_color_by_menu = False
        color_options = ['']

        scatter_data['_color'] = scatter_data['_black']
        scatter_data['_selection_color'] = scatter_data['_orange']
    
    else:
        show_color_by_menu = 'color_by' not in hide_widgets

        if isinstance(color_by, str):
            color_options = ['', color_by]
        else:
            color_options = [''] + list(color_by)

        scatter_data['_color'] = scatter_data[color_options[-1]]
        scatter_data['_selection_color'] = scatter_data[color_options[-1]]
    
    if label_by is None:
        label_by = df.index.name

    if isinstance(label_by, str):
        show_label_by_menu = False
        label_options = [label_by]
    else:
        show_label_by_menu = True
        label_options = list(label_by)
    
    scatter_data['_label'] = scatter_data[label_options[0]]

    if isinstance(marker_size, numbers.Number):
        size_widget_type = 'slider'
    else:
        if isinstance(marker_size, str):
            size_options = ['', marker_size]
        else:
            size_options = [''] + list(marker_size)
    
        size_widget_type = 'menu'
        marker_size = '_size'
        scatter_data[marker_size] = scatter_data[size_options[1]]

        scatter_data['_uniform_size'] = [6]*len(scatter_data[marker_size])

    scatter_source = bokeh.models.ColumnDataSource(data=scatter_data,
                                                   name='scatter_source',
                                                  )
    if len(initial_indices) > 0:
        scatter_source.selected = bokeh.models.Selection(indices=initial_indices)

    scatter = fig.scatter('x',
                          'y',
                          source=scatter_source,
                          size=marker_size,
                          fill_color='_color',
                          fill_alpha=initial_alpha,
                          line_color=None,
                          selection_color='_selection_color',
                          selection_alpha=0.9,
                          nonselection_color='_color',
                          nonselection_alpha=nonselection_alpha,
                          name='scatter',
                         )
    
    if log_scale:
        nonzero = df[numerical_cols][df[numerical_cols] > 0]

        overall_max = nonzero.max().max()
        overall_min = nonzero.min().min()
        
        initial = (overall_min * 0.1, overall_max * 10)
        bounds = (overall_min * 0.001, overall_max * 1000)
    
        def log(x):
            return np.log(x) / np.log(log_scale)

        bins = {}
        for name in numerical_cols:
            if identical_bins:
                left = overall_min * 0.9
                right = overall_max / 0.9
            else:
                name_min = nonzero[name].min()
                name_max = nonzero[name].max()
                left = name_min * 0.9
                right = name_max / 0.9

            bins[name] = list(np.logspace(log(left), log(right), num_bins, base=log_scale))

    else:
        overall_max = df.max(numeric_only=True).max()
        overall_min = df.min(numeric_only=True).min()
        overall_buffer = (overall_max - overall_min) * 0.05
        
        extent = overall_max - overall_min
        overhang = extent * 0.05
        max_overhang = extent * 0.5

        initial = (overall_min - overhang, overall_max + overhang)
        bounds = (overall_min - max_overhang, overall_max + max_overhang)
        
        bins = {}
        for name in numerical_cols:
            if identical_bins:
                left = overall_min - overall_buffer
                right = overall_max + overall_buffer
            else:
                name_min = df[name].min()
                name_max = df[name].max()
                name_buffer = (name_max - name_min) * 0.05
                left = name_min - name_buffer
                right = name_max + name_buffer
            
            bins[name] = list(np.linspace(left, right, num_bins))

    if initial_data_lims is not None:
        initial = initial_data_lims

    if data_bounds is not None:
        bounds = data_bounds
        # Ensure that auto-data limits don't extend out of bounds.
        initial = (max(initial[0], bounds[0]), min(initial[1], bounds[1]))

    diagonals_visible = ('diagonal' in grid)

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
    ) 

    diagonal = fig.line(x=bounds, y=bounds, **line_kwargs)
    lines = [
        diagonal,
    ]

    for line in lines:
        line.visible = diagonals_visible
    
    axes_lines = [
        fig.line(x=bounds, y=0, name='axes_line', **line_kwargs),
        fig.line(x=0, y=bounds, name='axes_line', **line_kwargs),
    ]

    for line in axes_lines:
        line.visible = show_axes_lines or ('axes' in grid)
    
    if volcano:
        fig.y_range = bokeh.models.Range1d(-0.1, 8)
        fig.x_range = bokeh.models.Range1d(-1, 1)
    else:
        if zoom_to_initial_data:
            x_min, x_max = bins[x_name][0], bins[x_name][-1]
            y_min, y_max = bins[y_name][0], bins[y_name][-1]
        else:
            x_min, x_max = initial
            y_min, y_max = initial

        fig.y_range = bokeh.models.Range1d(y_min, y_max)
        fig.x_range = bokeh.models.Range1d(x_min, x_max)

    fig.x_range.name = 'x_range'
    fig.y_range.name = 'y_range'
    
    lower_bound, upper_bound = bounds
    range_kwargs = dict(lower_bound=lower_bound, upper_bound=upper_bound)
    
    fig.outline_line_color = 'black'

    scatter.selection_glyph.line_color = None
    scatter.nonselection_glyph.line_color = None

    # Make marginal histograms.

    histogram_data = {}
    histogram_data = {
        'zero': [0]*(num_bins - 1),
    }

    for name in numerical_cols:
        histogram_data.update({
            f'{name}_bins_left': bins[name][:-1],
            f'{name}_bins_right': bins[name][1:],
        })

    max_count = 0
    for name in numerical_cols:
        counts, _ = np.histogram(df[name].dropna(), bins=bins[name])
        max_count = max(max(counts), max_count)
        histogram_data[f'{name}_all'] = list(counts)

    if log_scale:
        axis_type = 'log'
    else:
        axis_type = 'linear'

    hist_figs = {
        'x': bokeh.plotting.figure(width=plot_width, height=100,
                                   x_range=fig.x_range,
                                   x_axis_type=axis_type,
                                   name='hists_x',
                                  ),
        'y': bokeh.plotting.figure(width=100, height=plot_height,
                                   y_range=fig.y_range,
                                   y_axis_type=axis_type,
                                   name='hists_y',
                                   min_border_bottom=60 + (41 if two_level_index else 0), # mysterious behavior
                                  ),
    }

    if title is None:
        title = ''

    hist_figs['x'].title.text = title

    for axis, name in [('x', x_name), ('y', y_name)]:
        for data_type in ['all', 'bins_left', 'bins_right']:
            axis_key = f'{axis}_{data_type}'
            name_key = f'{name}_{data_type}'
            histogram_data[axis_key] = histogram_data[name_key]

        initial_vals = df[name].iloc[initial_indices]
        initial_counts, _ = np.histogram(initial_vals.dropna(), bins[name])
    
        histogram_data[f'{axis}_selected'] = initial_counts

    histogram_source = bokeh.models.ColumnDataSource(data=histogram_data,
                                                     name='histogram_source',
                                                    )
    initial_hist_alpha = 0.1 if len(initial_indices) > 0 else 0.2
    quads = {}
    quads['x_all'] = hist_figs['x'].quad(left='x_bins_left',
                                         right='x_bins_right',
                                         bottom='zero',
                                         top='x_all',
                                         source=histogram_source,
                                         color='black',
                                         alpha=initial_hist_alpha,
                                         line_color=None,
                                         name='hist_x_all',
                                        )

    quads['x_selected'] = hist_figs['x'].quad(left='x_bins_left',
                                              right='x_bins_right',
                                              bottom='zero',
                                              top='x_selected',
                                              source=histogram_source,
                                              color='orange',
                                              alpha=0.8,
                                              line_color=None,
                                             )

    quads['y_all'] = hist_figs['y'].quad(top='y_bins_left',
                                         bottom='y_bins_right',
                                         left='zero',
                                         right='y_all',
                                         source=histogram_source,
                                         color='black',
                                         alpha=initial_hist_alpha,
                                         line_color=None,
                                         name='hist_y_all',
                                        )

    quads['y_selected'] = hist_figs['y'].quad(top='y_bins_left',
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

    if identical_bins:
        x_end = max_count
        y_end = max_count
    else:
        x_end = max(histogram_source.data['x_all'])
        y_end = max(histogram_source.data['y_all'])

    hist_figs['x'].y_range = bokeh.models.Range1d(name='hist_x_range',
                                                  start=0,
                                                  end=x_end,
                                                  bounds='auto',
                                                 )
    hist_figs['y'].x_range = bokeh.models.Range1d(name='hist_y_range',
                                                  start=0,
                                                  end=y_end,
                                                  bounds='auto',
                                                 )

    for hist_fig in hist_figs.values():
        hist_fig.outline_line_color = None
        hist_fig.axis.visible = False
        hist_fig.grid.visible = False
        hist_fig.min_border = 0
        hist_fig.toolbar_location = None

    # Configure tooltips that pop up when hovering over a point.
    
    hover = bokeh.models.HoverTool(renderers=[scatter])
    hover.tooltips = [
        (df.index.name, '@{{{0}}}'.format(df.index.name)),
    ]
    for key in hover_keys:
        hover.tooltips.append((key, '@{{{0}}}'.format(key)))
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
        elif 'gene' in col_name:
            formatter = bokeh.models.widgets.HTMLTemplateFormatter(template='<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=<%= value %>" target="_blank"><%= value %></a>')
            width = 100
        else:
            formatter = bokeh.models.StringFormatter()
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
    
    filtered_source = bokeh.models.ColumnDataSource(data=filtered_data, name='filtered_source')
    
    table = bokeh.models.widgets.DataTable(source=filtered_source,
                                           columns=columns,
                                           width=2 * size if heatmap else size,
                                           height=600,
                                           sortable=False,
                                           reorderable=False,
                                           name='table',
                                           index_position=None,
                                          )
    
    labels = bokeh.models.LabelSet(x='x',
                                   y='y',
                                   text='_label',
                                   level='glyph',
                                   x_offset=0,
                                   y_offset=2,
                                   source=filtered_source,
                                   text_font_size=f'{label_size}pt',
                                   name='labels',
                                  )
    fig.add_layout(labels)
    
    # Set up menus or heatmap to select columns from df to put on x- and y-axis.

    if heatmap:
        correlations = df[numerical_cols].corr()
        correlations_array = np.array(correlations)
        if correlations_array.min() > 0:
            v_min = 0
            v_max = 1
            v_min = correlations_array.min()
            v_max = correlations_array.max()
            c_map = matplotlib.cm.viridis
        else:
            v_min = -1
            v_max = 1
            c_map = matplotlib.cm.RdBu_r

        norm = matplotlib.colors.Normalize(vmin=v_min, vmax=v_max)
        def r_to_color(r):
            color = matplotlib.colors.rgb2hex(c_map(norm(r)))
            return color
        
        linkage = scipy.cluster.hierarchy.linkage(correlations)
        dendro = scipy.cluster.hierarchy.dendrogram(linkage,
                                                    no_plot=True,
                                                    labels=numerical_cols,
                                                   )
        orders = {
            'original': numerical_cols,
            'clustered': dendro['ivl'],
        }
        data = defaultdict(list)

        for y, row in enumerate(numerical_cols):
            for x, col in enumerate(numerical_cols):
                for order_key, order in orders.items():
                    row = order[y]
                    col = order[x]

                    r = correlations[row][col]
                    data[f'r_{order_key}'].append(r)
                    data[f'x_{order_key}'].append(x)
                    data[f'x_name_{order_key}'].append(col)
                    data[f'y_{order_key}'].append(y)
                    data[f'y_name_{order_key}'].append(row)
                    data[f'color_{order_key}'].append(r_to_color(r))

        if cluster:
            order_key = 'clustered'
        else:
            order_key = 'original'

        for k in ['r', 'x', 'x_name', 'y', 'y_name', 'color']:
            data[k] = data[f'{k}_{order_key}']

        heatmap_source = bokeh.models.ColumnDataSource(data, name='heatmap_source')
        num_exps = len(numerical_cols)
        heatmap_size = size - 100
        heatmap_fig = bokeh.plotting.figure(tools='tap',
                                            x_range=(-0.5, num_exps - 0.5),
                                            y_range=(num_exps - 0.5, -0.5),
                                            width=heatmap_size, height=heatmap_size,
                                            toolbar_location=None,
                                            name='heatmap_fig',
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
        for axis in [heatmap_fig.xaxis, heatmap_fig.yaxis]:
            axis.name = 'heatmap_axis'

        hover = bokeh.models.HoverTool()
        hover.tooltips = [
            ('X', '@x_name'),
            ('Y', '@y_name'),
            ('r', '@r'),
        ]
        heatmap_fig.add_tools(hover)

        def make_tick_formatter(order):
            return '''var dict = {dict}; return dict[tick].slice(0, 15);'''.format(dict=dict(enumerate(order)))
        
        for ax in [heatmap_fig.xaxis, heatmap_fig.yaxis]:
            ax.ticker = bokeh.models.FixedTicker(ticks=np.arange(num_exps))
            ax.formatter = bokeh.models.FuncTickFormatter(code=make_tick_formatter(orders[order_key]))
            ax.major_label_text_font_size = '8pt'

        heatmap_fig.xaxis.major_label_orientation = np.pi / 4

        # Turn off black lines on bottom and left.
        for axis in (heatmap_fig.xaxis, heatmap_fig.yaxis):
            axis.axis_line_color = None

        name_pairs = list(zip(heatmap_source.data['x_name'], heatmap_source.data['y_name']))
        initial_index = name_pairs.index((x_name, y_name))
        heatmap_source.selected.indices = [initial_index]

        heatmap_fig.min_border = 1
        
        dendro_fig = bokeh.plotting.figure(height=100, width=heatmap_size,
                                           x_range=heatmap_fig.x_range,
                                          )
        icoord = np.array(dendro['icoord'])
        # Want to covert the arbitrary scaling that icoord is given in so that
        # the min value maps to 0 and the max value maps to len(correlations) - 1
        interval = (icoord.max() - icoord.min()) / (len(correlations) - 1)
        xs = list((icoord - icoord.min()) / interval)
        dcoord = np.array(dendro['dcoord'])
        ys = list(dcoord)
        dendrogram_lines = dendro_fig.multi_line(xs, ys, color='black')

        if not cluster:
            dendrogram_lines.visible = False

        dendro_fig.outline_line_color = None
        dendro_fig.axis.visible = False
        dendro_fig.grid.visible = False
        dendro_fig.min_border = 0
        dendro_fig.toolbar_location = None
        dendro_fig.y_range = bokeh.models.Range1d(start=0, end=dcoord.max() * 1.05)
    
        # Button to toggle clustering.
        cluster_button = bokeh.models.widgets.Toggle(label='cluster heatmap',
                                                     width=50,
                                                     name='cluster_button',
                                                     active=cluster,
                                                    )

    elif two_level_index:
        if not isinstance(menu_width, tuple):
            menu_width = (menu_width, menu_width)

        menus = {}

        for level, width in zip((0, 1), menu_width):
            level_name = original_df.columns.names[level]
            if level_name is None:
                level_name = level

            if level_name in level_order:
                options = level_order[level_name]
            else:
                options = [v for v in original_df.columns.levels[level] if v != 'gene' and v != '']

            for names, axis in zip(initial_xy_names, ('x', 'y')):
                menus[axis, level] = bokeh.models.widgets.MultiSelect(title=f'{axis.upper()}-axis: {level_name}',
                                                                      options=options,
                                                                      value=[names[level]],
                                                                      size=min(10, len(options)),
                                                                      name=f'{axis}_{level}_menu',
                                                                      width=width,
                                                                     )

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
        

    if not heatmap:
        # Just a placeholder so that we can assume cluster_button exists.
        cluster_button = bokeh.models.widgets.Toggle(name='cluster_button')

    # Button to toggle labels.
    label_button = bokeh.models.widgets.Toggle(label='Label selected points',
                                               width=50,
                                               active=True,
                                               name='label_button',
                                              )
    
    # Button to zoom to current data limits.
    zoom_to_data_button = bokeh.models.widgets.Button(label='Zoom to data limits',
                                                      width=50,
                                                      name='zoom_button',
                                                     )
    # Menu to choose label source.
    label_menu = bokeh.models.widgets.Select(title='Label by:',
                                             options=label_options,
                                             value=label_options[0],
                                             name='label_menu',
                                            )

    # Menu to choose color source.
    color_menu = bokeh.models.widgets.MultiSelect(title='Color by:',
                                                  options=color_options,
                                                  value=[color_options[-1]],
                                                  name='color_menu',
                                                 )

    # Radio group to choose whether to draw a vertical/horizontal grid or
    # diagonal guide lines. 
    options = ['grid', 'diagonal', 'diagonal+axes', 'none']
    active = options.index(grid)
    grid_options = bokeh.models.widgets.RadioGroup(labels=options,
                                                   active=active,
                                                   name='grid_radio_buttons',
                                                  )

    search_input = bokeh.models.widgets.TextInput(title='Search:')

    case_sensitive = bokeh.models.widgets.CheckboxGroup(labels=['Case sensitive'],
                                                        active=[],
                                                        name='case_sensitive',
                                                       )

    # Menu to select a subset of points from a columns of bools.
    subset_options = [''] + bool_cols
    subset_menu = bokeh.models.widgets.Select(title='Select subset:',
                                              options=subset_options,
                                              value='',
                                             )

    # button to dump table to file.
    save_button = bokeh.models.widgets.Button(label='Save table to file',
                                              width=50,
                                              name='save_button',
                                             )
    
    # button to dump annotations to file.
    annotation_button = bokeh.models.widgets.Button(label='save annotations to file',
                                                    width=50,
                                                    name='annotation_button',
                                                   )
    
    if alpha_widget_type == 'slider':
        alpha_widget = bokeh.models.Slider(start=0.,
                                           end=1.,
                                           value=initial_alpha,
                                           step=.05,
                                           title='alpha',
                                           name='alpha',
                                          )
    elif alpha_widget_type == 'text':
        alpha_widget = bokeh.models.TextInput(title='alpha', name='alpha', value=str(initial_alpha))
    else:
        raise ValueError('{0} not a valid alpha_widget_type value'.format(alpha_widget_type))

    if size_widget_type == 'slider':
        size_widget = bokeh.models.Slider(start=1,
                                          end=20.,
                                          value=marker_size,
                                          step=1,
                                          title='marker size',
                                          name='marker_size',
                                         )
    else:
        size_widget = bokeh.models.widgets.Select(title='Size by:',
                                                  options=size_options,
                                                  value=size_options[1],
                                                  name='marker_size',
                                                 )
    
    rect_source = bokeh.models.ColumnDataSource(data={
        'x': [],
        'y': [],
        'width': [],
        'height': [],
    }, name='rect_source')

    rects = fig.rect('x', 'y', 'width', 'height',
                     source=rect_source,
                     fill_alpha=0,
                     line_color='black',
                     line_alpha=0.5,
                    )
    rect_annotater = bokeh.models.BoxEditTool(renderers=[rects])
    if 'annotation' not in hide_widgets:
        fig.add_tools(rect_annotater)
    
    cross_source = bokeh.models.ColumnDataSource(data={
        'x': [],
        'y': [],
    }, name='cross_source')

    crosses = fig.cross('x', 'y', size=25,
                        line_color='black',
                        source=cross_source,
                       )
    cross_annotater = bokeh.models.PointDrawTool(renderers=[crosses])
    if 'annotation' not in hide_widgets:
        fig.add_tools(cross_annotater)

    fig.min_border = 1

    widgets = [
        label_button,
        zoom_to_data_button,
        cluster_button,
        label_menu,
        color_menu,
        alpha_widget,
        size_widget,
        grid_options,
        search_input,
        case_sensitive,
        subset_menu,
        save_button,
        annotation_button,
    ]

    if 'table' in hide_widgets:
        hide_widgets.add('save_button')

    if 'search' in hide_widgets:
        hide_widgets.add('case_sensitive')

    if 'annotation' in hide_widgets:
        hide_widgets.add('annotation_button')

    if not heatmap:
        hide_widgets.add('cluster_button')

    if not show_color_by_menu:
        hide_widgets.add('color_menu')

    if not show_label_by_menu:
        hide_widgets.add('label_menu')
    
    if len(subset_options) == 1:
        hide_widgets.add('subset_menu')

    widgets = [w for w in widgets if w.name not in hide_widgets]

    toolbar = bokeh.models.ToolbarBox(toolbar=fig.toolbar, toolbar_location='right')

    row = bokeh.layouts.row
    col = bokeh.layouts.column
    spacer = bokeh.layouts.Spacer

    columns = [
        col(hist_figs['x'], fig),
        col(spacer(height=100), hist_figs['y']),
        # Note: update to bokeh 2.4.1 changed behavior of ToolbarBox height
        # such that without this phantom 1px spacer, it thinks it can
        # only be 100px high and collapses most buttons into an overflow.
        col(spacer(height=100), row(spacer(width=1, height=800), toolbar)),
    ]

    if heatmap:
        heatmap_column = col(dendro_fig, heatmap_fig)
        columns.append(heatmap_column)
        widget_box = col(*widgets)
        final_column = col(spacer(height=min_border),
                           widget_box,
                          )
    else:
        if not two_level_index:
            widgets = [x_menu, y_menu] + widgets
            widget_box = col(*widgets)
            final_column = col(spacer(height=min_border),
                               widget_box,
                              )
        else:
            x_row = row(menus['x', 0], menus['x', 1])
            y_row = row(menus['y', 0], menus['y', 1])
            widget_box = col(*widgets)
            final_column = col(spacer(height=min_border),
                               x_row, y_row,
                               widget_box,
                              )

    columns.append(final_column)
        
    rows = [
        row(*columns),
    ]
    if 'table' not in hide_widgets:
        rows.append(table)

    full_layout = col(*rows)

    #rows = [
    #    hist_figs['x'],
    #    row(fig, hist_figs['y'], toolbar),
    #]

    #full_layout = col(*rows)


    # Build callbacks last so that all args are available.

    # Prevent ranges from going outside bounds.
    range_callback = build_js_callback(__file__, 'scatter_range',
                                       args=dict(
                                           x_range=fig.x_range,
                                           y_range=fig.y_range,
                                       ),
                                       format_kwargs=range_kwargs,
                                      )
    for range_ in [fig.x_range, fig.y_range]:
        for prop in ['start', 'end']:
            range_.js_on_change(prop, range_callback)

    # Filter the table when selection changes.
    selection_callback = build_js_callback(__file__, 'scatter_selection',
                                           args=dict(
                                               search_input=search_input,
                                               hist_x_all=quads['x_all'],
                                               hist_y_all=quads['y_all'],
                                               scatter_source=scatter_source,
                                               filtered_source=filtered_source,
                                               histogram_source=histogram_source,
                                               table=table,
                                               subset_menu=subset_menu
                                           ),
                                          )
    scatter_source.selected.js_on_change('indices', selection_callback)

    # Update displayed data when menu or heatmap selection changes.
    # Note: uses selection_callback, so must be after it.
    if heatmap:
        heatmap_callback = build_js_callback(__file__, 'scatter_heatmap',
                                             args=dict(
                                                 x_axis=fig.xaxis,
                                                 y_axis=fig.yaxis,
                                                 scatter_source=scatter_source,
                                                 filtered_source=filtered_source,
                                                 histogram_source=histogram_source,
                                                 heatmap_source=heatmap_source,
                                             ),
                                            )
        heatmap_source.selected.js_on_change('indices', heatmap_callback)

        # Cluster button
        format_kwargs = {
            # The callback expects code surrounded by double quotes.
            'clustered_formatter': '"{}"'.format(make_tick_formatter(orders['clustered'])),
            'original_formatter': '"{}"'.format(make_tick_formatter(orders['original'])),
        }
        cluster_callback = build_js_callback(__file__, 'scatter_cluster_button',
                                            format_kwargs=format_kwargs,
                                            args=dict(
                                                dendrogram_lines=dendrogram_lines,
                                                heatmap_source=heatmap_source,
                                                heatmap_x_axis=heatmap_fig.xaxis,
                                                heatmap_y_axis=heatmap_fig.yaxis,
                                                x_axis=fig.xaxis,
                                                y_axis=fig.yaxis,
                                                fig=fig,
                                            ),
                                            )
        cluster_button.js_on_click(cluster_callback)

    elif two_level_index:
        menu_callback = build_js_callback(__file__, 'scatter_two_level_menu',
                                          args=dict(
                                              scatter_source=scatter_source,
                                              filtered_source=filtered_source,
                                              histogram_source=histogram_source,
                                              x_menu_0=menus['x', 0],
                                              x_menu_1=menus['x', 1],
                                              y_menu_0=menus['y', 0],
                                              y_menu_1=menus['y', 1],
                                              x_axis_title=x_title,
                                              x_axis_subtitle=x_subtitle,
                                              y_axis_title=y_title,
                                              y_axis_subtitle=y_subtitle,

                                          ),
                                         )
        for axis in ('x', 'y'):
            for level in (0, 1):
                menus[axis, level].js_on_change('value', menu_callback)
    else:
        menu_callback = build_js_callback(__file__, 'scatter_menu',
                                          args=dict(
                                              scatter_source=scatter_source,
                                              filtered_source=filtered_source,
                                              histogram_source=histogram_source,
                                              x_menu=x_menu,
                                              y_menu=y_menu,
                                              x_axis=fig.xaxis,
                                              y_axis=fig.yaxis,
                                          ),
                                         )
        for menu in [x_menu, y_menu]:
            menu.js_on_change('value', menu_callback)

    format_kwargs = dict(log_scale=bool_to_js(log_scale),
                         identical_bins=bool_to_js(identical_bins),
                        )
    zoom_callback = build_js_callback(__file__, 'scatter_zoom_to_data',
                                      args=dict(
                                          scatter_source=scatter_source,
                                          histogram_source=histogram_source,
                                          x_range=fig.x_range,
                                          y_range=fig.y_range,
                                          hist_x_range=hist_figs['x'].y_range,
                                          hist_y_range=hist_figs['y'].x_range,
                                      ),
                                      format_kwargs=format_kwargs,
                                     )
    zoom_to_data_button.js_on_click(zoom_callback)

    # Toggle labels on points.
    label_button_callback = build_js_callback(__file__, 'scatter_label_button',
                                              args=dict(
                                                  labels=labels,
                                                  label_button=label_button,
                                              ),
                                             )
    label_button.js_on_click(label_button_callback)

    # Set point alpha.
    alpha_callback = build_js_callback(__file__, 'scatter_alpha',
                                       args=dict(
                                           scatter=scatter,
                                       ) 
                                      )
    alpha_widget.js_on_change('value', alpha_callback)


    # Update selection on text search or case sensitive toggle.
    columns_to_search = [c for c in object_cols if c not in color_options]
    search_callback = build_js_callback(__file__, 'scatter_search',
                                        args=dict(
                                            scatter_source=scatter_source,
                                            search_input=search_input,
                                            case_sensitive=case_sensitive,
                                        ),
                                        format_kwargs=dict(column_names=str(columns_to_search)),
                                       )
    search_input.js_on_change('value', search_callback)
    case_sensitive.js_on_change('active', search_callback)

    grid_callback = build_js_callback(__file__, 'scatter_grid',
                                      args=dict(
                                          fig_grid=fig.grid,
                                          fig_diagonal=diagonal,
                                          fig_axes_line=axes_lines,
                                      ),
                                     )
    grid_options.js_on_change('active', grid_callback)

    subset_callback = build_js_callback(__file__, 'scatter_subset_menu',
                                        format_kwargs=dict(subset_indices=str(subset_indices)),
                                        args=dict(
                                            scatter_source=scatter_source,
                                        )
                                       )
    subset_menu.js_on_change('value', subset_callback)

    # Marker size.
    if size_widget_type == 'slider':
        size_widget.js_link('value', scatter.glyph, 'size')
    else:
        size_callback = build_js_callback(__file__, 'scatter_size_menu',
                                          args=dict(
                                              scatter_source=scatter_source,
                                          ),
                                         )
        size_widget.js_on_change('value', size_callback)

    # Save button.
    save_callback = build_js_callback(__file__, 'scatter_save_button',
                                      format_kwargs=dict(column_names=str(table_col_names)),
                                      args=dict(
                                          filtered_source=filtered_source,
                                      ),
                                     )
    save_button.js_on_click(save_callback)

    # Label menu.
    label_callback = build_js_callback(__file__, 'scatter_label',
                                       args=dict(
                                           scatter_source=scatter_source,
                                           filtered_source=filtered_source,
                                       ),
                                      )
    label_menu.js_on_change('value', label_callback)

    # Color menu.
    color_callback = build_js_callback(__file__, 'scatter_color',
                                       args=dict(
                                           scatter_source=scatter_source,
                                       ),
                                      )
    color_menu.js_on_change('value', color_callback)

    # Save annotations.
    annotation_callback = build_js_callback(__file__, 'scatter_annotation_button',
                                            args=dict(
                                                rect_source=rect_source,
                                            ),
                                            format_kwargs=dict(
                                                title=title,
                                            ),
                                           )
    annotation_button.js_on_click(annotation_callback)
    
    if save_as is not None:
        bokeh.io.save(full_layout)
    else:
        bokeh.io.show(full_layout)

    if return_layout:
        return full_layout

def hex_to_CSS(hex_string, alpha=1.):
    ''' Converts an RGB hex value and optional alpha value to a CSS-format RGBA string. '''
    rgb = matplotlib.colors.colorConverter.to_rgb(hex_string)
    rgb = [int(v * 255) for v in rgb]
    CSS = 'rgba({1}, {2}, {3}, {0})'.format(alpha, *rgb)
    return CSS

def parallel_coordinates(df=None,
                         link_axes=True,
                         log_scale=False,
                         save_as=None,
                         initial_limits=None,
                         y_axis_label='',
                         min_start='-Infinity',
                         max_end='Infinity',
                         table_width=1200,
                         width=1200,
                         height=800,
                         label_by=None,
                         extra_injections=None,
                        ):
    ''' Makes an interactive parallel coordinates plot using d3. Call without
    any arguments for an example using data from Jan et al. Science 2014.
    Uses the parallel-coordinates library (github.com/syntagmatic/parallel-coordinates)
    Slickgrid integration taken from github.com/syntagmatic/parallel-coordinates/blob/master/examples/slickgrid.html
    Highlight-on-hover taken from bl.ocks.org/mostaphaRoudsari/b4e090bb50146d88aec4.

    Args:
        df: A pandas DataFrame with columns containing numerical data to plot. 
            Any text columns will be searchable through the 'Search:' field.
            A column named 'color' can be provided to color lines using any
            format understood by matplotlib.colors.to_rgb.
            If df is None, loads example data from Jan et al. Science 2014.

        link_axes: If True, all axes have the same range.

        log_scale: If link_axes=True, whether or not to use a (base 10) log scale.

        initial_limits: If link_axes=True, the intial start and end points of 
            axes limits.

        save_as: A filename (relative to the current directory, and probably 
            ending in .html) in which to save a standalone HTML file containing
            the figure. Opens a download dialog to save the file locally as well.
    '''
    if df is None:
        example_fn = Path(__file__).parent /  'jan_ratios.csv'
        df = pd.read_csv(example_fn, index_col='systematic_name')
        red = '#e41a1c'
        blue = '#0acbee'
        color = pd.Series('black', index=df.index)
        color[df['secretome']] = blue
        color[df['mitop2']] = red
        df['color'] = color

        log_scale = True

    if extra_injections is None:
        extra_injections = {}
    
    # Drop NaNs and copy before changing
    df = df.dropna().copy()
    
    # Collapse multiindex if present
    df.columns = [' '.join(map(str, n)) if isinstance(n, tuple) else str(n) for n in df.columns]

    # Collapsing multindex may mangle color column name.
    color_columns = [c for c in df.columns if c.startswith('color')]

    if len(color_columns) != 1:
        color_series = pd.Series('black', index=df.index)
    else:
        color_column = color_columns[0]
        color_series = df[color_column]
        df.drop(columns=color_column, inplace=True)

    def make_color_string(c):
        array = map(int, np.array(matplotlib.colors.to_rgb(c)) * 255)
        return 'rgba({0}, {1}, {2}, '.format(*array)
    df['_color'] = color_series.map(make_color_string)

    template_fn = Path(__file__).parent / 'template_inline.html'
    html_template = open(template_fn).read()

    encoded_data = base64.b64encode(df.to_csv().encode())
    URI = "'data:text/plain;base64,{0}'".format(encoded_data.decode())

    # Has to be checked in this order since bool subclasses int.
    if isinstance(log_scale, bool):
        log_scale = bool_to_js(log_scale)
    elif isinstance(log_scale, int):
        log_scale = str(log_scale)
    else:
        raise ValueError('log_scale should be an int or a bool', log_scale)

    # template_inline.html contains strings with format /*INJECT:KEY*/ at which data
    # corresponding to values of keys in injections dictionary below will be 
    # substituted into the html, replacing the entire comment-flanked string.
    injections = {
        'encoded_data': URI,
        'link_axes': bool_to_js(link_axes),
        'log_scale': log_scale,
        'initial_limits': str(list(initial_limits)) if initial_limits is not None else 'false',
        'y_axis_label': y_axis_label,
        'min_start': str(min_start),
        'max_end': str(max_end),
        'table_width': str(table_width),
        'width': str(width),
        'height': str(height),
        'label_by': str(label_by),
    }

    injections.update(extra_injections)

    def match_to_injection(match):
        # To accomodate optional injections, return the whole match if not registered.
        return injections.get(match.group(1), match.group(0))

    template_with_data = re.sub("\/\*INJECT:(.*)\*\/", match_to_injection, html_template)
    # Could probably find both /* and <-- in the same regex...
    template_with_data = re.sub("<!--*INJECT:(.*)-->", match_to_injection, template_with_data)
    
    if save_as is not None:
        with open(save_as, 'w') as fh:
            fh.write(template_with_data)
        output = None
    else:
        output = IPython.display.HTML(template_with_data)

    return output

def heatmap(df, cmap, vmin, vmax, height=500):
    normed = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)(df)
    rgba_float = cmap(normed)
    rgba_int = (rgba_float * 255).astype(np.uint8)[::-1] # flip row order for plotting

    # cryptic numpy manipulations convert from an NxNx4 array of uint8's
    # into an NxN array of uint32's.
    img = np.ascontiguousarray(rgba_int).view(dtype=np.uint32).reshape(df.shape)
    
    num_rows, num_cols = df.shape
    
    width = int(np.floor(height * num_cols / num_rows))
    fig = bokeh.plotting.figure(height=height, width=width, active_scroll='wheel_zoom')
    
    lower_bound = -0.5
    x_upper_bound = lower_bound + num_cols
    
    lower_bound = -0.5
    y_upper_bound = lower_bound + num_rows

    fig.image_rgba(image=[img], x=lower_bound, y=lower_bound, dw=num_cols, dh=num_rows)

    for axis in [fig.yaxis, fig.xaxis]:
        axis.ticker = np.arange(len(df))
        axis.major_label_text_font_size = '0pt'
        
    col_labels = [str(c) for c in df.columns.values]
    row_labels = [str(r) for r in df.index.values]

    fig.xaxis.major_label_overrides = {i: str(c) for i, c in enumerate(col_labels)}
    fig.yaxis.major_label_overrides = {i: l for i, l in enumerate(row_labels[::-1])}

    fig.xaxis.major_label_orientation = 'vertical'

    fig.x_range = bokeh.models.Range1d(lower_bound, x_upper_bound, name='x_range', tags=[x_upper_bound - lower_bound])
    fig.y_range = bokeh.models.Range1d(lower_bound, y_upper_bound, name='y_range')

    fig.grid.visible = False

    fig.axis.major_tick_line_color = None
    fig.axis.major_tick_out = 0
    fig.axis.major_tick_in = 0
    
    quad_source = bokeh.models.ColumnDataSource(data={
        'left': [],
        'right': [],
        'bottom': [],
        'top': [],
    }, name='quad_source',
    )
    fig.quad('left', 'right', 'top', 'bottom',
             fill_alpha=0,
             line_color='black',
             source=quad_source,
            )
    
    search_input = bokeh.models.TextInput(title='Search:', name='search')

    label_fig_size = 200
    label_figs = {
        'left': bokeh.plotting.figure(width=label_fig_size, height=height, x_range=(-10, 10), y_range=fig.y_range),
        'right': bokeh.plotting.figure(width=label_fig_size, height=height, x_range=(-10, 10), y_range=fig.y_range),
        'top': bokeh.plotting.figure(height=label_fig_size, width=width, y_range=(-10, 10), x_range=fig.x_range),
        'bottom': bokeh.plotting.figure(height=label_fig_size, width=width, y_range=(-10, 10), x_range=fig.x_range),
    }

    fig.min_border = 0
    fig.outline_line_color = None
    fig.axis.axis_line_color = None

    for label_fig in label_figs.values():
        label_fig.outline_line_color = None
        label_fig.axis.visible = False
        label_fig.grid.visible = False
        label_fig.min_border = 0
        label_fig.toolbar_location = None
    
    row_label_source = bokeh.models.ColumnDataSource(data={
        'text': row_labels,
        'smallest': [0]*len(row_labels),
        'biggest': [label_fig_size]*len(row_labels),
        'ascending': np.arange(len(row_labels)),
        'descending': np.arange(len(row_labels))[::-1],
        'alpha': [1]*len(row_labels),
    }, name='row_label_source')
    
    col_label_source = bokeh.models.ColumnDataSource(data={
        'text': col_labels,
        'smallest': [0]*len(col_labels),
        'biggest': [label_fig_size]*len(col_labels),
        'ascending': np.arange(len(col_labels)),
        'descending': np.arange(len(col_labels))[::-1],
        'alpha': [1]*len(col_labels),
    }, name='col_label_source')

    common_kwargs = {
        'text': 'text',
        'text_alpha': 'alpha',
        'text_font_size': f'{3000 // (max(x_upper_bound, y_upper_bound) - lower_bound)}%',
        'text_baseline': 'middle',
        'name': 'labels',
    }
    row_kwargs = {
        'y': 'descending',
        'x_units': 'screen',
        'y_units': 'data',
        'source': row_label_source, 
        **common_kwargs,
    }
    col_kwargs = {
        'x': 'ascending',
        'x_units': 'data',
        'y_units': 'screen',
        'source': col_label_source, 
        **common_kwargs,
    }

    labels = {
        'left': bokeh.models.LabelSet(x='biggest',
                                      text_align='right',
                                      x_offset=-5,
                                      **row_kwargs,
                                     ),
        'right': bokeh.models.LabelSet(x='smallest',
                                       text_align='left',
                                       x_offset=5,
                                       **row_kwargs,
                                      ),
        'top': bokeh.models.LabelSet(y='smallest',
                                     text_align='left',
                                     angle=np.pi / 2,
                                     y_offset=5,
                                     **col_kwargs,
                                    ),
        'bottom': bokeh.models.LabelSet(y='biggest',
                                        text_align='right',
                                        angle=np.pi / 2,
                                        y_offset=-5,
                                        **col_kwargs,
                                       ),
    }

    for which in labels:
        label_figs[which].add_layout(labels[which])
    
    fig.toolbar_location = None
    toolbar = bokeh.models.ToolbarBox(toolbar=fig.toolbar)

    # Build callbacks last so that all args are available.

    range_kwargs = dict(lower_bound=lower_bound, upper_bound=x_upper_bound)

    # Prevent ranges from going outside bounds, and resize fonts on zoom.
    range_callback = build_js_callback(__file__, 'heatmap_range',
                                       args=dict(
                                           labels=labels,
                                       ),
                                       format_kwargs=range_kwargs,
                                      )
    for range_ in [fig.x_range, fig.y_range]:
        for prop in ['start', 'end']:
            range_.js_on_change(prop, range_callback)

    search_format_kwargs = dict(
        row_labels=row_labels,
        col_labels=col_labels,
        lower_bound=lower_bound,
        x_upper_bound=x_upper_bound,
        y_upper_bound=y_upper_bound,
    )
    search_callback = build_js_callback(__file__, 'heatmap_search',
                                        args=dict(
                                            quad_source=quad_source,
                                            row_label_source=row_label_source,
                                            col_label_source=col_label_source,
                                        ),
                                        format_kwargs=search_format_kwargs,
                                       )
    search_input.js_on_change('value', search_callback)

    rows = [
        [bokeh.layouts.Spacer(width=label_fig_size), label_figs['top'], bokeh.layouts.Spacer(width=label_fig_size)],
        [label_figs['left'], fig, label_figs['right'], toolbar, search_input],
        [bokeh.layouts.Spacer(width=label_fig_size), label_figs['bottom']],
    ]

    layout = bokeh.layouts.column([bokeh.layouts.row(row) for row in rows])

    bokeh.io.show(layout)
