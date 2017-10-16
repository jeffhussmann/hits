import bokeh
from bokeh.models.tickers import SingleIntervalTicker
from bokeh.models.annotations import Legend, LegendItem
from bokeh.core.properties import List, Instance
import copy
import json
import os.path
import numpy as np
import Sequencing.genetic_code as genetic_code
from itertools import cycle
from collections import defaultdict
from external_coffeescript import build_callback
from colors_list import colors_list

# In 0.12.5 and 0.12.6, custom models have to be defined BEFORE output_notebook
# is called. See https://github.com/codypiersall/bokeh-custom-model.
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

bokeh.io.reset_output()
bokeh.io.output_notebook()

def codon(enrichments=None,
          group_by='codon',
          groupings=None,
          colors=None,
          y_max=5,
          x_lims=(-25, 25),
          unselected_alpha=0.2,
          initial_menu_selection=None,
          initial_top_group_selections=None,
          initial_sub_group_selections=None,
          intial_resolution='codon',
         ):
    ''' An interactive plot of metacodon enrichment profiles using bokeh. Call
    without any arguments for an example using data from Jan et al. Science
    2014.

    Args:

        enrichments: A multi-level dictionary of enrichment values to plot, laid
            out like:
            
            enrichments = {
                'codon': {
                    'xs': [list of codon offset values],
                    'experiment_1': {
                        'TTT': [list of codon-resolution enrichments],
                        'TTC': ...,
                        ...,
                    },
                    'experiment_2': {...},
                },
                'nucleotide': {
                    'xs': [list of nucleotide offset values],
                    'experiment_1': {
                        'TTT': [list of nucleotide-resolution enrichments],
                        'TTC': ...,
                        ...,
                    },
                    'experiment_2': {...},
                },
            }

            See example_metacodon_enrichments.json in the same directory as this
            file for an example. (If enrichments == None, loads this file, which
            contains processed data from Jan et al. Science 2014, as an example.
            Other arguments are overwritten to highlight interesting features.)
        
        group_by: If 'codon', plot enrichments around one codon for all
            experiments. If 'experiment', plot enrichments around all codons for
            one experiment.

        groupings: If not None, a dictionary of groups of experiments/codons to
            list together for convenient selection.
            If None and group_by == 'experiment', experiment names will be
            grouped by name.split(':')[0].
            If None and group_by == 'codon', codons will be grouped by amino
            acid.
    '''
    if initial_top_group_selections is None:
        initial_top_group_selections = []
    if initial_sub_group_selections is None:
        initial_sub_group_selections = []
    
    if enrichments is None:
        # Load example data
        dirname = os.path.dirname(__file__)
        fn = os.path.join(dirname, 'example_metacodon_enrichments.json')
        with open(fn) as fh:
            enrichments = json.load(fh)

        # Prepare example groupings.
        exp_names = set(enrichments['codon'].keys())
        exp_names.remove('xs')
        group_names = ['-CHX', '+CHX_2min', '+CHX_7min']
        groupings = {}
        for group_name in group_names:
            groupings[group_name] = [n for n in exp_names if group_name in n]

        # Overrule other arguments to highlight interesting features in data.
        x_lims = (-18, 81)
        y_max = 3.2
        group_by = 'codon'
        initial_menu_selection = 'CGA'
        initial_sub_group_selections = [
            'BirA_+CHX_2minBiotin_input',
            'BirAmVenusUbc6_+CHX_7minBiotin_input',
            'sec63mVenusBirA_-CHX_1minBiotin_input',
        ]

    # enrichments has experiments as top level keys and codons as lower level
    # keys. Building ColumnDataSource's below assumes a dictionary with
    # checkbox_names as the top keys, so if grouping by codon, don't need
    # to do anything. If grouping by experiment, need to invert the order of the
    # dictionary.

    # Copy enrichments since xs will be popped.
    enrichments = copy.deepcopy(enrichments)

    xs = {
        'codon': enrichments['codon'].pop('xs'),
        'nucleotide': enrichments['nucleotide'].pop('xs'),
    }
    
    exp_names = sorted(enrichments['codon'])

    if group_by == 'experiment':
        if groupings is None:
            groupings = {aa: cs for aa, cs in genetic_code.full_back_table.items() if aa != '*'}

        menu_options = exp_names
        checkbox_names = genetic_code.non_stop_codons
        
        inverted = {}
        for resolution in enrichments:
            inverted[resolution] = {}
            for checkbox_name in checkbox_names:
                inverted[resolution][checkbox_name] = {}
                for menu_name in menu_options:
                    inverted[resolution][checkbox_name][menu_name] = enrichments[resolution][menu_name][checkbox_name]

        enrichments = inverted
    
    elif group_by == 'codon':
        menu_options = []
        for codon in genetic_code.non_stop_codons:
            aa = genetic_code.forward_table[codon]
            codon_aa = '{0} ({1})'.format(codon, aa)
            menu_options.append(codon_aa)
            for exp_name in exp_names:
                for resolution in enrichments:
                    ys = enrichments[resolution][exp_name].pop(codon)
                    enrichments[resolution][exp_name][codon_aa] = ys

        checkbox_names = sorted(enrichments['codon'])
        
        if groupings is None:
            groupings = defaultdict(list)
            for checkbox_name in sorted(checkbox_names):
                group_name = checkbox_name.split(':')[0]
                groupings[group_name].append(checkbox_name)

    if colors is None:
        colors = dict(zip(checkbox_names, cycle(colors_list)))

    if initial_menu_selection is None:
        initial_menu_selection = menu_options[0]
    if initial_menu_selection not in menu_options:
        raise ValueError('{0} not in {1}'.format(initial_menu_selection, menu_options))

    # Select all members of any initially selected top_group. 
    initial_sub_group_selections = set(initial_sub_group_selections)
    for top_name, sub_names in sorted(groupings.items()):
        if top_name in initial_top_group_selections:
            initial_sub_group_selections.update(sub_names)
    
    # Build ColumnDataSource's from enrichments.
    sources = {}
    for key in ['plotted', 'codon', 'nucleotide']:
        if key == 'plotted':
            resolution = intial_resolution
        else:
            resolution = key
        
        sources[key] = {}
        for checkbox_name in sorted(enrichments[resolution]):
            source = bokeh.models.ColumnDataSource(enrichments[resolution][checkbox_name])
            source.data['x'] = xs[resolution]
            source.data['y'] = source.data[initial_menu_selection]
            source.data['name'] = [checkbox_name] * len(xs[resolution])
            source.name = 'source_{0}_{1}'.format(checkbox_name, key)
            sources[key][checkbox_name] = source
   

    # Set up the actual plot.
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
    fig.toolbar.logo = None

    fig.grid.grid_line_alpha = 0.4

    fig.y_range = bokeh.models.Range1d(0, y_max)
    fig.y_range.name = 'y_range'
    fig.x_range = bokeh.models.Range1d(*x_lims)
    fig.x_range.name = 'x_range'

    range_callback = build_callback('metacodon_range')
    fig.y_range.callback = range_callback
    fig.x_range.callback = range_callback

    fig.xaxis.axis_label = 'Offset ({0}s)'.format(intial_resolution)
    fig.yaxis.axis_label = 'Mean relative enrichment'

    fig.xaxis.name = 'x_axis'
    fig.yaxis.name = 'y_axis'

    fig.xaxis.axis_label_text_font_style = 'normal'
    fig.yaxis.axis_label_text_font_style = 'normal'
    
    fig.xaxis[0].ticker = bokeh.models.tickers.SingleIntervalTicker(interval=3, num_minor_ticks=3)

    legend_items = []
    initial_legend_items = []
    lines = []
    for checkbox_name, source in sources['plotted'].items():
        if checkbox_name in initial_sub_group_selections:
            color = colors[checkbox_name]
            line_width = 2
            line_alpha = 0.95
            circle_alpha = 0.9
        else:
            color ='black'
            line_width = 1
            circle_alpha = 0
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
                            size=4.5,
                            fill_alpha=circle_alpha,
                            line_alpha=0,
                            hover_alpha=0.95,
                            hover_color=colors[checkbox_name],
                           )
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

    source_callback = build_callback('metacodon_selection')
    for source in sources['plotted'].values():
        source.callback = source_callback

    hover = bokeh.models.HoverTool(line_policy='interp',
                                   renderers=lines,
                                  )
    hover.tooltips = [('name', '@name')]
    fig.add_tools(hover)

    # Draw horizontal and vertical lines.
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

    menu_title = 'Codon:' if group_by == 'codon' else 'Experiment:'
    menu = bokeh.models.widgets.MultiSelect(options=menu_options,
                                            value=[initial_menu_selection],
                                            size=min(30, len(menu_options)),
                                            title=menu_title,
                                           )
    menu.callback = build_callback('metacodon_menu')

    sub_group_callback = build_callback('metacodon_sub_group',
                                        format_kwargs=dict(color_unselected='false'),
                                       )

    top_group_callback = build_callback('metacodon_top_group')

    top_groups = []
    sub_groups = []

    for top_name, sub_names in sorted(groupings.items()):
        width = int(75 + max(len(l) for l in sub_names) * 6.5)

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

    highest_level_chooser = bokeh.models.widgets.RadioGroup(labels=['codon resolution', 'nucleotide resolution'],
                                                            active=0 if intial_resolution == 'codon' else 1,
                                                            name='highest_level_chooser',
                                                           )

    injection_sources = []
    for resolution in ['codon', 'nucleotide']:
        injection_sources.extend(sources[resolution].values())
    injection = {'ensure_no_collision_{0}'.format(i): v for i, v in enumerate(injection_sources)}

    highest_level_chooser.callback = build_callback('metacodon_resolution',
                                                    args=injection,
                                                   )

    clear_selection = bokeh.models.widgets.Button(label='Clear selection')
    clear_selection.callback = build_callback('metacodon_clear_selection')
    
    alpha_slider = bokeh.models.Slider(start=0.,
                                       end=1.,
                                       value=unselected_alpha,
                                       step=.05,
                                       title='unselected alpha',
                                      )
    alpha_slider.callback = build_callback('lengths_unselected_alpha')

    widgets = [
        menu,
        highest_level_chooser,
        alpha_slider,
        clear_selection,
    ]

    grid = [
        top_groups,
        sub_groups,
        [fig, bokeh.layouts.widgetbox(widgets)],
    ]

    bokeh.io.show(bokeh.layouts.layout(grid))

def gene(enrichments=None,
         groupings=None,
         y_max=5,
         x_lims=(-20, 150),
         unselected_alpha=0.2,
         initial_assignment=None,
         initial_top_group_selections=None,
         initial_sub_group_selections=None,
         landmark_pairs=None,
         group_order=None,
         colors=None,
        ):
    ''' An interactive plot of metagene enrichment profiles using bokeh. Call
    without any arguments for an example using data from Jan et al. Science
    2014.

    Args:

        enrichments: A multi-level dictionary of enrichment values to plot, laid
            out like:
            
            enrichments = {
                'codon': {
                    'start_codon': {
                        'xs': [list of codon offset values relative to start codon],
                        'experiment_1': [list of enrichment values],
                        'experiment_2': [...],
                        ...,
                    },
                    'stop_codon': {...}
                },
                'nucleotide': {...}
            }

            See example_metagene_enrichments.json in the same directory as this
            file for an example. (If enrichments == None, loads this file, which
            contains processed data from Jan et al. Science 2014, as an 
            example.)
        
        groupings: If not None, a dictionary of groups of experiments to list
            together for convenient selection.
            If None, experiment names will be grouped by name.split(':')[0].
    '''
    if initial_top_group_selections is None:
        initial_top_group_selections = []
    if initial_sub_group_selections is None:
        initial_sub_group_selections = []

    if landmark_pairs is None:
        landmark_pairs = [('start_codon', 'stop_codon')]
    
    if enrichments is None:
        # Load example data
        dirname = os.path.dirname(__file__)
        fn = os.path.join(dirname, 'example_metagene_enrichments.json')
        with open(fn) as fh:
            enrichments = json.load(fh)

        # Prepare example groupings.
        exp_names = set(enrichments['codon']['start_codon'].keys())
        exp_names.remove('xs')
        group_names = ['-CHX', '+CHX_2min', '+CHX_7min']
        groupings = {}
        for group_name in group_names:
            groupings[group_name] = [n for n in exp_names if group_name in n]

        # Overrule other arguments to highlight interesting features in data.
        initial_resolution = 'codon'
        y_max = 18
        initial_sub_group_selections = [
            'BirA_+CHX_2minBiotin_input',
            'BirAmVenusUbc6_+CHX_7minBiotin_input',
            'sec63mVenusBirA_-CHX_1minBiotin_input',
        ]
    
    enrichments = copy.deepcopy(enrichments)

    assignments = sorted(enrichments)
    assignment = assignments[0]
    if initial_assignment is None:
        initial_assignment = assignment
    landmarks = sorted(enrichments[assignment])

    all_xs = {}
    for assignment in assignments:
        all_xs[assignment] = {}
        for landmark in landmarks:
            all_xs[assignment][landmark] = enrichments[assignment][landmark].pop('xs')
    
    exp_names = sorted(enrichments[assignment][landmarks[0]])

    if groupings is None:
        groupings = defaultdict(list)
        for exp_name in exp_names:
            group_name = exp_name.split(':')[0]
            groupings[group_name].append(exp_name)

    initial_sub_group_selections = set(initial_sub_group_selections)
    for top_name, sub_names in sorted(groupings.items()):
        if top_name in initial_top_group_selections:
            initial_sub_group_selections.update(sub_names)

    sources = {}

    landmark_pair_names = []
    for left, right in landmark_pairs:
        landmark_pair_name = ':'.join((left, right))
        landmark_pair_names.append(landmark_pair_name)
        for assignment in assignments:
            sources[assignment, landmark_pair_name] = {}

            for exp_name in exp_names:
                source = bokeh.models.ColumnDataSource()
                
                for side, landmark in (('left', left), ('right', right)):
                    source.data['xs_{0}'.format(side)] = all_xs[assignment][landmark]
                    source.data['ys_{0}'.format(side)] = enrichments[assignment][landmark][exp_name]
                
                source.data['name'] = [exp_name] * len(source.data['xs_left'])
            
                source.name = 'source_{0}_{1}_{2}'.format(exp_name, assignment, landmark_pair_name)
                sources[assignment, landmark_pair_name][exp_name] = source

    sources['plotted'] = {}
    initial_sources = sources[initial_assignment, landmark_pair_names[0]]
    for exp_name in exp_names:
        initial_data = initial_sources[exp_name].data
        source = bokeh.models.ColumnDataSource(data=initial_data)
        source.name = 'source_{0}_plotted'.format(exp_name)
        sources['plotted'][exp_name] = source

    initial_lims = {
        'left': x_lims,
        'right': -1 * np.array(x_lims[::-1]),
    }

    x_ranges = {}
    x_range_callback = build_callback('metagene_x_range')
    for side in ('left', 'right'):
        x_range = bokeh.models.Range1d(*initial_lims[side])
        x_range.name = 'x_range_{0}'.format(side)
        x_range.callback = x_range_callback
        x_ranges[side] = x_range

    y_range = bokeh.models.Range1d(0, y_max)
    y_range.name = 'y_range'
    y_range.callback = build_callback('metagene_y_range')

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
    for side in ('left', 'right'):
        figs[side] = bokeh.plotting.figure(plot_width=800,
                                           plot_height=600,
                                           x_range=x_ranges[side],
                                           y_range=y_range,
                                           y_axis_location=side,
                                           tools=tools,
                                           active_scroll='wheel_zoom',
                                          )

    lines = {
        'left': [],
        'right': [],
    }

    legend_items = []
    initial_legend_items = []

    if colors is None:
        colors = dict(zip(exp_names, cycle(colors_list)))

    for exp_name in exp_names:
        if exp_name in initial_sub_group_selections:
            color = colors[exp_name]
            line_width = 2
            line_alpha = 0.95
            circle_alpha = 0.9
        else:
            color = 'black'
            line_width = 1
            circle_alpha = 0
            if len(initial_sub_group_selections) > 0:
                line_alpha = unselected_alpha
            else:
                line_alpha = 0.6

        for side in ('left', 'right'):
            line = figs[side].line(x='xs_{0}'.format(side),
                                   y='ys_{0}'.format(side),
                                   source=sources['plotted'][exp_name],
                                   color=color,
                                   nonselection_line_color=colors[exp_name],
                                   nonselection_line_alpha=unselected_alpha,
                                   hover_alpha=1.0,
                                   hover_color=colors[exp_name],
                                   line_width=line_width,
                                   line_alpha=line_alpha,
                                  )
            line.hover_glyph.line_width = 3
            line.name = 'line_{0}'.format(exp_name)
            lines[side].append(line)

            circle = figs[side].circle(x='xs_{0}'.format(side),
                                       y='ys_{0}'.format(side),
                                       source=sources['plotted'][exp_name],
                                       size=4,
                                       color=colors[exp_name],
                                       fill_alpha=circle_alpha,
                                       line_alpha=circle_alpha,
                                       hover_alpha=1.0,
                                       hover_color=colors[exp_name],
                                      )
            #circle.hover_glyph.visible = True
            circle.name = 'circle_{0}'.format(exp_name)

            if side == 'right':
                legend_item = LegendItem(label=exp_name, renderers=[line])
                legend_items.append(legend_item)
                if exp_name in initial_sub_group_selections:
                    initial_legend_items.append(legend_item)

    legend = ToggleLegend(name='legend',
                          items=initial_legend_items,
                          all_items=legend_items,
                         )
    figs['right'].add_layout(legend)
    figs['right'].legend.location = 'top_left'
    figs['right'].legend.background_fill_alpha = 0.5
    
    for side, landmark in zip(('left', 'right'), landmark_pairs[0]):
        fig = figs[side]
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
        
        remove_underscore = ' '.join(landmark.split('_'))
        fig.xaxis.axis_label = 'Offset from {0}'.format(remove_underscore)
        fig.xaxis.axis_label_text_font_style = 'normal'
        fig.xaxis.axis_label_text_font_size = '16pt'
        
        fig.grid.grid_line_alpha = 0.4

    figs['left'].yaxis.axis_label = 'Mean relative enrichment'
    figs['left'].yaxis.axis_label_text_font_style = 'normal'
    figs['left'].yaxis.axis_label_text_font_size = '16pt'
    
    source_callback = build_callback('metacodon_selection')
    for source in sources['plotted'].values():
        source.callback = source_callback

    for side in ['left', 'right']:
        hover = bokeh.models.HoverTool(line_policy='interp',
                                       renderers=lines[side],
                                      )
        hover.tooltips = [('name', '@name')]
        figs[side].add_tools(hover)

    injection_sources = []
    for assignment in assignments:
        for landmark_pair_name in landmark_pair_names:
            to_inject = sources[assignment, landmark_pair_name].values()
            injection_sources.extend(to_inject)

    injection = {'ensure_no_collision_{0}'.format(i): v for i, v in enumerate(injection_sources)}

    assignment_menu = bokeh.models.widgets.Select(title='Assignment:',
                                                  options=assignments,
                                                  value=initial_assignment,
                                                 )
    assignment_menu.name = 'assignment_menu'
    
    assignment_menu.callback = build_callback('metacodon_assignment',
                                              args=injection,
                                             )
        
    landmark_menu = bokeh.models.widgets.Select(title='Landmarks:',
                                                options=landmark_pair_names,
                                                value=landmark_pair_names[0],
                                               )
    landmark_menu.name = 'landmark_menu'
    
    landmark_menu.callback = build_callback('metacodon_assignment',
                                            args=injection,
                                           )
    
    clear_selection = bokeh.models.widgets.Button(label='Clear selection')
    clear_selection.callback = build_callback('metacodon_clear_selection')
    
    sub_group_callback = build_callback('metacodon_sub_group',
                                        format_kwargs=dict(color_unselected='false'),
                                       )

    top_group_callback = build_callback('metacodon_top_group')

    top_groups = []
    sub_groups = []
    if group_order is None:
        group_order = sorted(groupings)

    for top_name in group_order:
        sub_names = sorted(groupings[top_name])
        width = int(75 + max(len(l) for l in sub_names) * 6.5)
        
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
    
    alpha_slider = bokeh.models.Slider(start=0.,
                                       end=1.,
                                       value=unselected_alpha,
                                       step=.05,
                                       title='unselected alpha',
                                      )
    alpha_slider.callback = build_callback('lengths_unselected_alpha')

    plots = bokeh.layouts.gridplot([[figs['left'], figs['right']]])
    plots.children[0].logo = None

    grid = [
        top_groups,
        sub_groups,
        [plots,
         bokeh.layouts.widgetbox([assignment_menu,
                                  landmark_menu,
                                  alpha_slider,
                                  clear_selection,
                                 ]),
        ],
    ]
    bokeh.io.show(bokeh.layouts.layout(grid))

def lengths(raw_counts,
            group_by='experiment',
            groupings=None,
            initial_top_group_selections=None,
            initial_sub_group_selections=None,
            types=None,
            max_length=1000,
            x_max=500,
            colors=None,
           ):
    ys = {
        'raw_counts': raw_counts,
        'by_type': {},
        'by_exp': {},
    }

    for exp_name in ys['raw_counts']:
        total_counts = sum(ys['raw_counts'][exp_name]['insert_lengths'])

        ys['by_type'][exp_name] = {}
        ys['by_exp'][exp_name] = {}

        for type_name in raw_counts[exp_name]:
            raw = raw_counts[exp_name][type_name]
            ys['by_type'][exp_name][type_name] = np.true_divide(raw, sum(raw))
            ys['by_exp'][exp_name][type_name] = np.true_divide(raw, total_counts)

    first_exp_name = ys['raw_counts'].keys()[0]
    
    if types is None:
        types = sorted(ys['raw_counts'][first_exp_name])


    # ys has experiments as top level keys and types as lower level keys.
    # Need a dictionary with checkbox names as the top keys, so if grouping by
    # type, don't need to do anything. If grouping by experiment, need to invert
    # the order.
    # 2017_03_06: this comment doesnt seem to make sense.

    if group_by == 'experiment':
        menu_options = sorted(ys['raw_counts'])
        checkbox_names = types
        
        rekeyed_ys = {}
        for normalization in ys:
            rekeyed_ys[normalization] = {}
            for checkbox_name in checkbox_names:
                rekeyed_ys[normalization][checkbox_name] = {}
                for menu_name in menu_options:
                    rekeyed_ys[normalization][checkbox_name][menu_name] = ys[normalization][menu_name][checkbox_name]

        ys = rekeyed_ys

    elif group_by == 'type':
        menu_options = types
        checkbox_names = sorted(ys['raw_counts'])

    if initial_sub_group_selections is None:
        initial_sub_group_selections = []

    if initial_top_group_selections is None:
        initial_top_group_selections = []

    if groupings is None:
        groupings = {n: [n] for n in checkbox_names}

    if colors is None:
        colors = dict(zip(checkbox_names, cycle(colors_list)))

    initial_menu_selection = None
    if initial_menu_selection is None:
        initial_menu_selection = menu_options[0]

    if initial_menu_selection not in menu_options:
        raise ValueError('{0} not in {1}'.format(initial_menu_selection, menu_options))

    sources = {}
    for key in ['plotted'] + ys.keys():
        if key == 'plotted':
            normalization = 'raw_counts'
        else:
            normalization = key

        sources[key] = {}
        for checkbox_name in checkbox_names:
            data = {}
            for menu_name in ys[normalization][checkbox_name]:
                full_list = list(ys[normalization][checkbox_name][menu_name])
                if len(full_list) > max_length:
                    full_list[max_length] = sum(full_list[max_length:])
                elif len(full_list < max_length):
                    full_list.extend([0]*(max_length + 1 - len(full_list)))

                data[menu_name] = np.array(full_list[:max_length + 1])
            source = bokeh.models.ColumnDataSource(data)

            xs = np.arange(max_length + 1)
            source.data['x'] = xs

            source.data['y'] = source.data[initial_menu_selection]

            source.data['name'] = [checkbox_name] * len(xs)
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

    fig.y_range.name = 'y_range'
    fig.x_range.name = 'x_range'
    fig.x_range.end = x_max

    range_callback = build_callback('lengths_range')
    fig.y_range.callback = range_callback
    fig.x_range.callback = range_callback

    fig.yaxis.axis_label = 'Number of reads'
    fig.yaxis.axis_label_text_font_style = 'normal'
    fig.yaxis.axis_label_text_font_size = '16pt'

    fig.xaxis.name = 'x_axis'
    fig.yaxis.name = 'y_axis'

    fig.xaxis.axis_label = 'Insert length'
    fig.xaxis.axis_label_text_font_style = 'normal'
    fig.xaxis.axis_label_text_font_size = '16pt'

    fig.yaxis.axis_label_text_font_style = 'normal'

    legend_items = []
    initial_legend_items = []
    lines = []
    for checkbox_name, source in sources['plotted'].items():
        if checkbox_name in initial_sub_group_selections:
            color = colors[checkbox_name]
            line_width = 2
            line_alpha = 0.95
            circle_alpha = 0.9
        else:
            color = colors[checkbox_name]
            line_width = 1
            circle_alpha = 0
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
                            size=1,
                            fill_alpha=circle_alpha,
                            line_alpha=circle_alpha,
                            hover_alpha=0.95,
                            hover_color=colors[checkbox_name],
                           )
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
    
    source_callback = build_callback('metacodon_selection')
    for source in sources['plotted'].values():
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
    menu.callback = build_callback('metacodon_menu')
    
    sub_group_callback = build_callback('metacodon_sub_group',
                                        format_kwargs=dict(color_unselected='true'),
                                       )
    top_group_callback = build_callback('metacodon_top_group')

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
    
    alpha_slider = bokeh.models.Slider(start=0.,
                                       end=1.,
                                       value=0.5,
                                       step=.05,
                                       title='alpha',
                                      )
    alpha_slider.callback = build_callback('lengths_unselected_alpha')
    
    highest_level_chooser = bokeh.models.widgets.Select(options=ys.keys(),
                                                        value='raw_counts',
                                                       )

    injection_sources = []
    for key in ys:
        injection_sources.extend(sources[key].values())
    # Note: use throwaway names instead of source names to ensure no illegal
    # characters in names.
    injection = {'ensure_no_collision_{0}'.format(i): v for i, v in enumerate(injection_sources)}

    highest_level_chooser.callback = build_callback('lengths_highest_level',
                                                    args=injection,
                                                   )

    grid = [
        top_groups,
        sub_groups,
        [fig, bokeh.layouts.widgetbox([menu, highest_level_chooser, alpha_slider])],
    ]
    bokeh.io.show(bokeh.layouts.layout(grid))
