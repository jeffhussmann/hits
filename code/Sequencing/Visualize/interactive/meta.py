import bokeh
import positions
from . import external_coffeescript

colors_list = [ # from http://godsnotwheregodsnot.blogspot.ru/2012/09/color-distribution-methodology.html
    "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
    "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
    "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
    "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
    "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
    "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
    "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
    "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
    "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
    "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
]

def codon(xs, ys, colors, groupings,
          toggle_resolution=True,
          y_max=5,
          x_max=25,
          unselected_alpha=0.2,
         ):
    sources = {}

    highest_level_keys = sorted(ys)

    for key in ['plotted'] + highest_level_keys:
        if key == 'plotted':
            resolution = highest_level_keys[0]
        else:
            resolution = key
        
        sources[key] = {}
        for checkbox_name in sorted(ys[resolution]):
            menu_names = sorted(ys[resolution][checkbox_name])
            source = bokeh.models.ColumnDataSource(ys[resolution][checkbox_name])
            source.data['x'] = xs[resolution]
            source.data['y'] = source.data[menu_names[0]]
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

    fig = bokeh.plotting.figure(plot_width=1600, plot_height=800,
                                tools=tools, active_scroll='wheel_zoom',
                               )

    fig.y_range = bokeh.models.Range1d(0, y_max)
    fig.x_range = bokeh.models.Range1d(-x_max, x_max)
    fig.x_range.name = 'x_range'

    range_callback =  external_coffeescript('metacodon_range',
                                            args=dict(fig=fig),
                                           )
    fig.y_range.callback = range_callback
    fig.x_range.callback = range_callback

    legend_items = []
    lines = []
    for checkbox_name, source in sources['plotted'].items():
        line = fig.line(x='x',
                        y='y',
                        color='black',
                        source=source,
                        line_width=1,
                        line_alpha=0.6,
                        line_join='round',
                        hover_alpha=1.0,
                        hover_color=colors[checkbox_name],
                        legend=checkbox_name,
                       )
        line.hover_glyph.line_width = 4
        line.name = 'line_{0}'.format(checkbox_name)
        lines.append(line)
        
        circle = fig.circle(x='x',
                            y='y',
                            color='black',
                            source=source,
                            size=4,
                            fill_alpha=0.9,
                            line_alpha=0.9,
                            visible=False,
                            hover_alpha=1.0,
                            hover_color=colors[checkbox_name],
                           )
        circle.hover_glyph.visible = True
        circle.name = 'circle_{0}'.format(checkbox_name)
    
        legend_items.append((checkbox_name, [line]))
        
    fig.legend.name = 'legend'
    fig.legend.items = []
        
    invisible_legend = bokeh.models.Legend(items=legend_items, name='invisible_legend')
    
    source_callback = external_coffeescript('metacodon_selection',
                                            args=dict(invisible_legend=invisible_legend),
                                           )
    for source in sources['plotted'].values():
        source.callback = source_callback

    hover = bokeh.models.HoverTool(line_policy='interp',
                                   renderers=lines,
                                  )
    hover.tooltips = [('name', '@name')]
    fig.add_tools(hover)

    zero = bokeh.models.annotations.Span(location=0,
                                         dimension='height',
                                         line_color='black',
                                         line_alpha=0.8,
                                         line_dash='dashed',
                                        )
    fig.renderers.append(zero)

    options = sorted(ys[highest_level_keys[0]].values()[0].keys())
    menu = bokeh.models.widgets.Select(options=options, value=options[0])
    menu.callback = external_coffeescript('metacodon_menu')

    sub_group_callback = external_coffeescript('metacodon_sub_group',
                                               format_args=dict(colors_dict=colors,
                                                                unselected_alpha=unselected_alpha
                                                               ),
                                               args=dict(invisible_legend=invisible_legend),
                                              )

    top_group_callback = external_coffeescript('metacodon_top_group',
                                               args=dict(invisible_legend=invisible_legend),
                                              )

    top_groups = []
    sub_groups = []
    width = 75 + max(len(l) for top_name in groupings for l in groupings[top_name]) * 5
    for top_name, sub_names in sorted(groupings.items()):
        top = bokeh.models.widgets.CheckboxGroup(labels=[top_name],
                                                 active=[],
                                                 width=width,
                                                 name='top_{0}'.format(top_name),
                                                 callback=top_group_callback,
                                                )
        top_groups.append(top)
        sub = bokeh.models.widgets.CheckboxGroup(labels=sub_names,
                                                 active=[],
                                                 width=width,
                                                 callback=sub_group_callback,
                                                 name='sub_{0}'.format(top_name),
                                                )
        sub_groups.append(sub)

    if toggle_resolution:
        highest_level_chooser = bokeh.models.widgets.RadioGroup(labels=['codon resolution', 'nucleotide resolution'], active=0)
        callback_name = 'metacodon_resolution'
    else:
        highest_level_chooser = bokeh.models.widgets.Select(options=highest_level_keys,
                                                            value=highest_level_keys[0],
                                                           )
        callback_name = 'metacodon_highest_level'
    
    injection_sources = []
    for key in highest_level_keys:
        injection_sources.extend(sources[key].values())
    injection = {'ensure_no_collision_{0}'.format(i): v for i, v in enumerate(injection_sources)}

    highest_level_chooser.callback = external_coffeescript(callback_name,
                                                           args=dict(fig=fig, **injection),
                                                          )

    clear_selection = bokeh.models.widgets.Button(label='Clear selection')
    clear_selection.callback = external_coffeescript('metacodon_clear_selection',
                                                     args=dict(invisible_legend=invisible_legend),
                                                    )

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
        ):
    exp_names = sorted(densities['codon'])
    sources = {}

    highest_level_keys = sorted(densities)
    
    max_before = 90
    max_after = 250
    xs_dict = {
        'codon': {'start_codon': np.arange(-max_before, max_after),
                  'stop_codon': np.arange(-max_after, max_before),
                 },
        'nucleotide': {'start_codon': np.arange(-3 * max_before, 3 * max_after),
                       'stop_codon': np.arange(-3 * max_after, 3 * max_before),
                      },
    }

    for key in ['plotted'] + highest_level_keys:
        if key == 'plotted':
            resolution = highest_level_keys[0]
        else:
            resolution = key
        
        sources[key] = {}
        for exp_name in exp_names:
            source = bokeh.models.ColumnDataSource()
            for landmark in ['start_codon', 'stop_codon']:
                xs = xs_dict[resolution][landmark]
                density_type = positions.MetageneDensityType('three_prime',
                                                            landmark,
                                                            'all',
                                                            'none_nearby',
                                                            0.1,
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
        figs[key] = bokeh.plotting.figure(plot_width=600,
                                          plot_height=600,
                                          x_range=x_ranges[key],
                                          y_range=y_range,
                                          y_axis_location=y_axis_location,
                                          tools=tools,
                                          active_scroll='wheel_zoom',
                                         )

    lines = {'start_codon': [], 'stop_codon': []}
    legend_items = []

    colors = dict(zip(exp_names, bokeh.palettes.brewer['Dark2'][8]))
    for exp_name in exp_names:
        for landmark in ('start_codon', 'stop_codon'):
            if landmark == 'stop_codon':
                legend_kwarg = {'legend': exp_name}
            else:
                legend_kwarg = {}

            line = figs[landmark].line(x='xs_{0}'.format(landmark),
                                       y='ys_{0}'.format(landmark),
                                       source=sources['plotted'][exp_name],
                                       color=colors[exp_name],
                                       hover_alpha=1.0,
                                       hover_color=colors[exp_name],
                                       line_width=1.5,
                                       line_alpha=0.6,
                                       **legend_kwarg)
            line.hover_glyph.line_width = 4
            line.name = 'line_{0}'.format(exp_name)
            lines[landmark].append(line)
            
            circle = figs[landmark].circle(x='xs_{0}'.format(landmark),
                                           y='ys_{0}'.format(landmark),
                                           source=sources['plotted'][exp_name],
                                           size=4,
                                           color=colors[exp_name],
                                           fill_alpha=0.9,
                                           line_alpha=0.9,
                                           visible=False,
                                           hover_alpha=1.0,
                                           hover_color=colors[exp_name],
                                          )
            circle.hover_glyph.visible = True
            circle.name = 'circle_{0}'.format(exp_name)

            if landmark == 'stop_codon':
                legend_items.append((exp_name, [line]))
    
    figs['stop_codon'].legend.name = 'legend'
    figs['stop_codon'].legend.items = []
    
    invisible_legend = bokeh.models.Legend(items=legend_items, name='invisible_legend')

    source_callback = external_coffeescript('metacodon_selection',
                                            args=dict(invisible_legend=invisible_legend),
                                           )
    for source in sources['plotted'].values():
        source.callback = source_callback

    for landmark in ['start_codon', 'stop_codon']:
        hover = bokeh.models.HoverTool(line_policy='interp',
                                       renderers=lines[landmark],
                                      )
        hover.tooltips = [('name', '@name')]
        figs[landmark].add_tools(hover)

    resolution = bokeh.models.widgets.RadioGroup(labels=['codon resolution', 'nucleotide resolution'], active=0)
    resolution.name = 'resolution'
    
    injection_sources = []
    for key in highest_level_keys:
        injection_sources.extend(sources[key].values())
    injection = {'ensure_no_collision_{0}'.format(i): v for i, v in enumerate(injection_sources)}

    resolution.callback = external_coffeescript('metacodon_resolution',
                                       args=dict(**injection),
                                      )
    
    sub_group_callback = external_coffeescript('metacodon_sub_group',
                                               format_args=dict(colors_dict=colors,
                                                                unselected_alpha=unselected_alpha
                                                               ),
                                               args=dict(invisible_legend=invisible_legend),
                                              )

    top_group_callback = external_coffeescript('metacodon_top_group',
                                               args=dict(invisible_legend=invisible_legend),
                                              )

    top_groups = []
    sub_groups = []
    width = 75 + max(len(l) for top_name in groupings for l in groupings[top_name]) * 6
    for top_name, sub_names in sorted(groupings.items()):
        top = bokeh.models.widgets.CheckboxGroup(labels=[top_name],
                                                 active=[],
                                                 width=width,
                                                 name='top_{0}'.format(top_name),
                                                 callback=top_group_callback,
                                                )
        top_groups.append(top)
        sub = bokeh.models.widgets.CheckboxGroup(labels=sub_names,
                                                 active=[],
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
