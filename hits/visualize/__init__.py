import colorsys
import io
import itertools
import functools
from collections import Counter

import scipy.stats
import seaborn as sns
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import PIL

from .. import utilities
from . import define_igv_colors

igv_colors = define_igv_colors.normalized_rgbs

blues_cdict = {
    'red':   ((0.0, 1.0, 1.0),
              (1.0, 0.0, 0.0)),
    'green': ((0.0, 1.0, 1.0),
              (1.0, 0.0, 0.0)),
    'blue':  ((0.0, 1.0, 1.0),
              (1.0, 1.0, 1.0)),
}

blues = matplotlib.colors.LinearSegmentedColormap('blues', blues_cdict, 10000)
blues.set_over('black')

reds_cdict = {
    'red':   ((0.0, 1.0, 1.0),
              (1.0, 1.0, 1.0)),
    'green': ((0.0, 1.0, 1.0),
              (1.0, 0.0, 0.0)),
    'blue':  ((0.0, 1.0, 1.0),
              (1.0, 0.0, 0.0)),
}

reds = matplotlib.colors.LinearSegmentedColormap('reds', reds_cdict, 10000)
reds.set_over('black')

greens_cdict = {
    'red':   ((0.0, 1.0, 1.0),
              (1.0, 0.0, 0.0)),
    'green': ((0.0, 1.0, 1.0),
              (1.0, 1.0, 1.0)),
    'blue':  ((0.0, 1.0, 1.0),
              (1.0, 0.0, 0.0)),
}

greens = matplotlib.colors.LinearSegmentedColormap('greens', greens_cdict, 10000)

def optional_ax(original_function):
    @functools.wraps(original_function)
    def possibly_new_ax(*args, **kwargs):
        ax_given = kwargs.get('ax')
        if not ax_given:
            fig, ax = plt.subplots(figsize=(12, 8))
            kwargs['ax'] = ax
        
        result = original_function(*args, **kwargs)
        
        figure_file_name = kwargs.get('save_as')
        if figure_file_name:
            assert not ax_given, 'Can\'t give ax and figure_file_name'
            fig.savefig(figure_file_name, bbox_inches='tight')
            plt.close(fig)
        
        return result
    
    return possibly_new_ax

def add_commas_to_ticks(ax, which='y'):
    def commas_formatter(x, pos):
        return f'{int(x):,}'

    tick_formatter = matplotlib.ticker.FuncFormatter(commas_formatter)
    if which == 'y':
        axis = ax.yaxis
    elif which == 'x':
        axis = ax.xaxis
    else:
        raise ValueError(which)

    axis.set_major_formatter(tick_formatter)

@optional_ax
def enhanced_scatter(xs, ys,
                     ax=None,
                     data=None,
                     colors=None,
                     do_fit=False,
                     show_correlation=False,
                     show_p_value=False,
                     show_beta=True,
                     hists_location=None,
                     hist_bins=100,
                     hist_range=None,
                     hist_alpha=0.2,
                     hist_size_ratio=0.1,
                     hist_colors=None,
                     hist_share_max=True,
                     marker_size=4,
                     text_size=14,
                     text_weight='normal',
                     text_location='lower right',
                     color_by_correlation=False,
                     fit_line_kwargs={'color': 'black',
                                      'alpha': 0.3,
                                     },
                     r_digits=2,
                     variance=False,
                     in_log_space=False,
                     big_fit_line=False,
                     colorbar=False,
                     remove_y_hist=False,
                     remove_x_hist=False,
                     label=None,
                     alpha=None,
                     **scatter_kwargs,
                    ):

    if data is not None:
        x_label = xs
        y_label = ys
        color_label = colors
        xs = data[xs]
        ys = data[ys]
        try:
            if colors in data:
                colors = data[colors]
        except TypeError:
            pass
    else:
        x_label = ''
        y_label = ''
        color_label = ''

    xs = np.asarray(xs)
    ys = np.asarray(ys)

    same_lists = np.allclose(xs, ys)

    if isinstance(colors, str) and colors == 'density':
        if not same_lists and len(xs) > 2:
            indices = np.arange(len(xs))
            np.random.shuffle(indices)
            random_indices = indices[:5000]

            sampled_points = np.vstack([xs[random_indices], ys[random_indices]])
            points = np.vstack([xs, ys])
            kernel = scipy.stats.gaussian_kde(sampled_points)
            colors = kernel(points)
        else:
            colors = 'black'

    elif colors is not None:
        if isinstance(colors, str):
            # colors is a single color name
            pass
        else:
            #colors = np.asarray(colors)
            pass
    
    if same_lists:
        do_fit = False

    kwargs = {
        's': marker_size,
        'label': label,
        'alpha': alpha,
        **scatter_kwargs,
    }

    kwargs.setdefault('linewidths', (0,))

    scatter = ax.scatter(xs, ys, c=colors, **kwargs)
    ax.set_xlabel(x_label, size=16)
    ax.set_ylabel(y_label, size=16)
    fig = ax.figure
    ax_position = ax.get_position()

    cax = None

    if not isinstance(colors, str):
        if colorbar:
            cax = fig.add_axes((ax_position.x0 + ax_position.width * 1.17,
                                ax_position.y1 - ax_position.height * 0.4,
                                ax_position.width * 0.03,
                                ax_position.height * 0.35,
                               ),
                              )
            fig.colorbar(scatter, cax=cax)
            cax.set_title(color_label)

    if isinstance(text_location, tuple):
        x, y = text_location
        horizontal_alignment = 'center'
        vertical_alignment = 'top'
        x_sign = 1
        y_sign = -1

    if 'left' in text_location:
        x = 0
        x_sign = 1
        horizontal_alignment = 'left'
    elif 'right' in text_location:
        x = 1
        x_sign = -1
        horizontal_alignment = 'right'

    if 'upper' in text_location:
        y = 1
        y_sign = -1
        vertical_alignment = 'top'
    elif 'lower' in text_location:
        y = 0
        y_sign = 1
        vertical_alignment = 'bottom'

    if text_location == 'above':
        x = 0.5
        y = 1
        x_sign = 1
        y_sign = 1.05
        vertical_alignment = 'bottom'
        horizontal_alignment = 'center'
        x_offset = 0
        y_offset = 0
    else:
        x_offset = 10
        y_offset = 0.5 * text_size

    text_kwargs = {
        'xy': (x, y),
        'xycoords': 'axes fraction',
        'textcoords': 'offset points',
        'horizontalalignment': horizontal_alignment,
        'verticalalignment': vertical_alignment,
        'fontsize': text_size,
        'family': 'serif',
        'weight': text_weight,
    }
    
    if do_fit:
        fit = np.polyfit(xs, ys, 1)
        beta, constant = fit
        fit_function = np.poly1d(fit)
        x_lims = ax.get_xlim()
        if big_fit_line:
            plot_lims = (min(xs) - (max(xs) - min(xs)), max(xs) + (max(xs) - min(xs)))
        else:
            plot_lims = x_lims
        ax.plot(plot_lims, fit_function(plot_lims), **fit_line_kwargs)
        ax.set_xlim(*x_lims)
        
        if show_beta == 'fit':
            ax.annotate(r'$y = {0:0.2f} x {2:s} {1:0.2f}$'.format(beta, abs(constant), '+' if constant > 0 else '-'),
                        xytext=(x_sign * x_offset, y_sign * y_offset * 4),
                        **text_kwargs)
        elif show_beta:
            ax.annotate(r'$\beta$ = {:0.2f}'.format(beta),
                        xytext=(x_sign * x_offset, y_sign * y_offset * 4),
                        **text_kwargs)

    if show_correlation:
        if in_log_space:
            # Values have been log'ed before being passed in, so exponentiate them
            # to recover correct correlation.
            r, p = scipy.stats.pearsonr(in_log_space**xs, in_log_space**ys)
        else:
            r, p = scipy.stats.pearsonr(xs, ys)

        if variance:
            r_part = 'r$^2$' + '= {:0.{digits}f}'.format(r**2, digits=r_digits)
        else:
            r_part = 'r = {:0.{digits}f}'.format(r, digits=r_digits)

        if show_p_value:
            p_part = ', p={:0.2e}'.format(p)
        else:
            p_part = ''

        text = r_part + p_part

        if color_by_correlation:
            text_kwargs['color'] = matplotlib.cm.seismic(0.5 * r + 0.5)

        ax.annotate(text,
                    xytext=(x_sign * x_offset, y_sign * y_offset),
                    **text_kwargs,
                   )

    max_x = None
    max_y = None

    if hists_location is not None:
        if hists_location == 'inside':
            bottom = ax_position.y0
            left = ax_position.x0
        elif hists_location == 'outside':
            bottom = ax_position.y1
            left = ax_position.x1
        
        common_kwargs = {
            'alpha': hist_alpha,
            'histtype': 'stepfilled',
        }

        if isinstance(hist_bins, int):
            hist_bins = {which: hist_bins for which in ['x', 'y']} 

        if hist_colors is None:
            hist_colors = {which: colors if isinstance(colors, str) else 'black' for which in ['x', 'y']}

        if hist_range is None:
            hist_range = {
                'x': (min(xs), max(xs)),
                'y': (min(ys), max(ys)),
            }

        ax_x = fig.add_axes((ax_position.x0, bottom, ax_position.width, ax_position.height * hist_size_ratio), sharex=ax)
        n_x, *rest = ax_x.hist(xs, range=hist_range['x'], bins=hist_bins['x'], color=hist_colors['x'], **common_kwargs)
        ax_x.axis('off')

        ax_y = fig.add_axes((left, ax_position.y0, ax_position.width * hist_size_ratio, ax_position.height), sharey=ax)
        n_y, *rest = ax_y.hist(ys, range=hist_range['y'], orientation='horizontal', bins=hist_bins['y'], color=hist_colors['y'], **common_kwargs)
        ax_y.axis('off')

        if hist_share_max:
            max_n = max(max(n_x), max(n_y))
            max_x = max_n
            max_y = max_n
        else:
            max_x = max(n_x)
            max_y = max(n_y)

        ax_x.set_ylim(0, max_x * 1.1 if max_x > 0 else 1)
        ax_y.set_xlim(0, max_y * 1.1 if max_y > 0 else 1)

        if remove_x_hist:
            fig.delaxes(ax_x)
        if remove_y_hist:
            fig.delaxes(ax_y)
    else:
        ax_x = None
        ax_y = None

    return fig, {'scatter': ax, 'hist_x': ax_x, 'hist_y': ax_y, 'colorbar': cax, 'hist_maxes': {'x': max_x, 'y': max_y}}

def draw_diagonal(ax, anti=False, color='black', **kwargs):
    if anti:
        xs, ys = [0, 1], [1, 0]
    else:
        xs, ys = [0, 1], [0, 1]
        if ax.get_xlim() != ax.get_ylim():
            print('warning: diagonal in non-equal axes')
    
    ax.plot(xs, ys,
            transform=ax.transAxes,
            color=color,
            **kwargs)

def draw_zeros_and_diagonal(ax, color='black', alpha=0.5, **kwargs):
    draw_diagonal(ax, color=color, alpha=alpha, **kwargs)
    ax.axhline(0, color=color, alpha=alpha, **kwargs)
    ax.axvline(0, color=color, alpha=alpha, **kwargs)

def label_scatter_plot(ax, xs, ys, labels,
                       data=None,
                       color=None,
                       to_label=slice(None),
                       vector='orthogonal',
                       initial_distance=5,
                       distance_increment=10,
                       arrow_alpha=0.2,
                       manual_ratios=None,
                       manual_alignments=None,
                       text_kwargs={'size': 10},
                       avoid=True,
                       avoid_axis_labels=False,
                       avoid_existing=False,
                       min_arrow_distance=10,
                       slope=1,
                      ):
    if data is not None:
        xs = data[xs]
        ys = data[ys]

        if color is not None:
            if color in data:
                color = data[color]
            else:
                color = [color]*len(xs)

        if isinstance(labels, list):
            pass
        elif labels in data:
            labels = data[labels]
        elif labels == data.index.name:
            labels = data.index
        else:
            raise IndexError

    def attempt_text(x, y, site, distance, vector, color):
        if vector == 'orthogonal':
            x_offset = np.sign(x * slope - y) * distance
            y_offset = -np.sign(x * slope - y) * distance
            ha = 'center'
            va = 'top' if y_offset < 0 else 'bottom'
        elif vector == 'upper left':
            x_offset = -distance
            y_offset = distance
            ha = 'center'
            va = 'bottom'
        elif vector == 'lower left':
            x_offset = -distance
            y_offset = -distance
            ha = 'center'
            va = 'top'
        elif vector == 'lower right':
            x_offset = distance
            y_offset = -distance
            ha = 'center'
            va = 'top'
        elif vector == 'upper right':
            x_offset = distance
            y_offset = distance
            ha = 'center'
            va = 'bottom'
        elif vector == 'radial':
            norm = np.linalg.norm([x, y])
            x_offset = x * distance / norm
            y_offset = y * distance / norm
            ha = 'center'
            va = 'top' if y_offset < 0 else 'bottom'
        elif vector == 'above':
            x_offset = 0
            y_offset = distance
            ha = 'center'
            va = 'bottom'
        elif vector == 'below':
            x_offset = 0
            y_offset = -distance
            ha = 'center'
            va = 'top'
        elif vector == 'sideways':
            x_offset = distance
            y_offset = 0
            ha, va = 'left', 'center'
        elif vector == 'right':
            x_offset = distance
            y_offset = 0
            ha, va = 'left', 'center'
        elif vector == 'left':
            x_offset = -distance
            y_offset = 0
            ha, va = 'right', 'center'
        elif vector == 'manual':
            x_ratio, y_ratio = manual_ratios
            ha, va = manual_alignments
            x_offset = distance * x_ratio
            y_offset = distance * y_ratio
        else:
            raise ValueError(vector)

        text = ax.annotate(site,
                           xy=(x, y),
                           xycoords=('data', 'data'),
                           xytext=(x_offset, y_offset),
                           textcoords='offset points',
                           ha=ha,
                           va=va,
                           color=color,
                           **text_kwargs,
                          )

        if avoid:
            ax.figure.canvas.draw()

        return text, text.get_window_extent(), (x, y, x_offset, y_offset)

    ax.figure.canvas.draw()

    starting_labels = []

    if avoid_axis_labels:
        to_add = [ax.xaxis.get_label(), ax.yaxis.get_label()] + ax.get_yticklabels() + ax.get_xticklabels()
        starting_labels.extend(to_add)

    if avoid_existing:
        to_add = [c for c in ax.get_children() if isinstance(c, matplotlib.text.Annotation)]
        starting_labels.extend(to_add)

    bboxes = []
    for label in starting_labels:
        try:
            bboxes.append(label.get_window_extent())
        except:
            pass

    if isinstance(vector, str):
        vector = np.array([vector]*len(xs))

    if color is None:
        color = 'black'

    if isinstance(color, str):
        color = np.array([color]*len(xs))

    tuples = zip(
        xs[to_label],
        ys[to_label],
        labels[to_label],
        vector[to_label],
        color[to_label],
    )

    for x, y, label, vec, color in tuples:
        distance = initial_distance
        text, bbox, coords = attempt_text(x, y, label, distance, vec, color)

        while avoid and any(bbox.overlaps(other_bbox) for other_bbox in bboxes):
            text.remove()
            distance += distance_increment
            text, bbox, coords = attempt_text(x, y, label, distance, vec, color)
            if distance >= distance_increment * 50:
                print(f'gave up on {label}')
                break
        
        if distance >= min_arrow_distance:
            x, y, x_offset, y_offset = coords
            ax.annotate('',
                        xy=(x, y),
                        xycoords=('data', 'data'),
                        xytext=(x_offset, y_offset),
                        textcoords=('offset points', 'offset points'),
                        arrowprops={'arrowstyle': '-',
                                    'alpha': arrow_alpha,
                                    'color': color,
                                   },
                       )

        bboxes.append(bbox)

def evenly_spaced_jet_colors(n):
    jet_colors = [matplotlib.cm.jet(x) for x in np.linspace(0, 1, n)]
    return jet_colors

def assign_colors_by_ranks(key_to_ranks, index=0):
    colors = evenly_spaced_jet_colors(len(key_to_ranks))
    key_to_color = {key: colors[key_to_ranks[key][index]] for key in key_to_ranks}
    return key_to_color

def color_labels(labels, name_to_color):
    for label in labels:
        color = name_to_color.get(label.get_text(), None)
        if color:
            label.set_color(color)

def apply_alpha(color, alpha, multiplicative=False):
    r, g, b, a =  matplotlib.colors.colorConverter.to_rgba(color)
    if multiplicative:
        a *= alpha
    else:
        a = alpha
    return (r, g, b, a)

def scale_darkness(color, scale):
    ''' Increases darkness of color by converting to hsl and dividing l value
    by scale. adapted from https://stackoverflow.com/a/49601444
    '''
    rgb = matplotlib.colors.to_rgb(color)
    h, l, s = colorsys.rgb_to_hls(*rgb)
    scaled_rgb = colorsys.hls_to_rgb(h, max(0, min(1, l / scale)), s)
    scaled_hex = matplotlib.colors.to_hex(scaled_rgb)
    
    return scaled_hex

def force_integer_ticks(axis):
    axis.set_major_locator(matplotlib.ticker.MaxNLocator(interger=True))

@optional_ax
def plot_counts(l, ax=None, log_scales=None, normalize=False, **kwargs):
    if log_scales is None:
        log_scales = set()

    counts = Counter(l)
    ys = utilities.counts_to_array(counts)

    first_nonzero = ys.nonzero()[0][0]
    xs = np.arange(first_nonzero, len(ys))
    ys = ys[first_nonzero:]
    if normalize:
        ys = ys / ys.sum()

    ax.plot(xs, ys, **kwargs)

    if 'x' in log_scales:
        ax.set_xscale('log')
    else:
        ax.set_xlim(-0.01 * len(ys), 1.01 * len(ys))

    if 'y' in log_scales:
        ax.set_yscale('log')
    else:
        ax.set_ylim(0)

def make_stacked_Image(figs, orientation='vertical', dpi=None):
    ims = []

    for fig in figs:
        with io.BytesIO() as buffer:
            fig.savefig(buffer, format='png', bbox_inches='tight', dpi=dpi)
            im = PIL.Image.open(buffer)
            im.load()
            ims.append(im)

        plt.close(fig)
        
    if not ims:
        return None

    if orientation == 'vertical':
        total_height = sum(im.height for im in ims)
        max_width = max(im.width for im in ims)

        stacked_im = PIL.Image.new('RGBA', size=(max_width, total_height), color='white')
        y_start = 0
        for im in ims:
            stacked_im.paste(im, (max_width - im.width, y_start))
            y_start += im.height

    else:
        max_height = max(im.height for im in ims)
        total_width = sum(im.width for im in ims)

        stacked_im = PIL.Image.new('RGBA', size=(total_width, max_height), color='white')
        x_start = 0
        for im in ims:
            stacked_im.paste(im, (x_start, max_height - im.height))
            x_start += im.width

    return stacked_im

def assign_categorical_colors(series, palette=None, sort=True):
    if palette is None:
        palette = sns.color_palette()

    values = series.unique()

    if sort:
        values = sorted(values)

    value_to_color = dict(zip(values, itertools.cycle(palette)))
    colors = series.map(value_to_color)

    return colors, value_to_color

def draw_categorical_legend(value_to_color,
                            ax,
                            font_size=12,
                            legend_location='upper left',
                            order=None,
                            manual_xy=None,
                            aliases=None,
                            title=None,
                           ):

    if legend_location == 'upper left':
        xy = (0.04, 0.96)
        va = 'top'
        ha  = 'left'
    elif legend_location == 'upper middle':
        xy = (0.5, 0.96)
        va = 'top'
        ha = 'center'
    elif legend_location == 'middle left':
        xy = (0.04, 0.5)
        va = 'center'
        ha = 'left'
    elif legend_location == 'upper right':
        xy = (0.96, 0.96)
        va = 'top'
        ha = 'right'
    elif legend_location == 'lower right':
        xy = (0.96, 0.4)
        va = 'bottom'
        ha = 'right'
    elif legend_location == 'lower left':
        xy = (0.06, 0.4)
        va = 'bottom'
        ha = 'left'
    elif legend_location == 'middle':
        xy = (0.5, 0.5)
        va = 'center'
        ha = 'center'

    if manual_xy is not None:
        xy = manual_xy

    if order is None:
        order = sorted(value_to_color)

    if aliases is None:
        aliases = {}

    if title != None:
        order = [f'{title}:'] + order

    for i, category in enumerate(order):
        if category == f'{title}:':
            color = 'black'
        else:
            color = value_to_color[category]
        ax.annotate(aliases.get(category, category),
                    xy=xy,
                    xycoords='axes fraction',
                    xytext=(0, -(font_size + 3) * i),
                    textcoords='offset points',
                    color=color,
                    va=va,
                    ha=ha,
                    size=font_size,
                   )

def add_y_axis_inset_zoom_ax(main_ax, inset_y_lim, location='right'):
    fig = main_ax.figure

    inverted_fig_tranform = fig.transFigure.inverted().transform    

    def draw_line(path, **kwargs):
        path_in_fig = [inverted_fig_tranform(ax.get_yaxis_transform().transform((x, y))) for x, y, ax in path]
        
        xs = [point[0] for point in path_in_fig]
        ys = [point[1] for point in path_in_fig]
        
        line = matplotlib.lines.Line2D(xs,
                                       ys, 
                                       transform=fig.transFigure,
                                       clip_on=False,
                                       color='black',
                                       solid_capstyle='butt',
                                       **kwargs,
                                      )
        fig.lines.append(line)
    
    main_brackset_sign = -1

    if location == 'right':
        ax_bounds = [1.2, 0, 1, 1]
        main_bracket_x = 1.04
        inset_bracket_x = -0.1
        inset_bracket_sign = 1

    else:
        ax_bounds = [0, -1.1, 1, 1]
        main_bracket_x = 1.05
        inset_bracket_x = 1.02
        inset_bracket_sign = -1

    inset_ax = main_ax.inset_axes(ax_bounds, sharex=main_ax)

    inset_ax.set_ylim(*inset_y_lim)

    bracket_width = 0.01
    bracket_offset = 0
    
    paths = [
        [
            (main_bracket_x, inset_y_lim[1], main_ax),
            (inset_bracket_x, inset_y_lim[1], inset_ax),
        ],
        [
            (main_bracket_x, inset_y_lim[0], main_ax),
            (inset_bracket_x, inset_y_lim[0], inset_ax),
        ],
    ]
    
    for path in paths:
        draw_line(path, alpha=0.5, linestyle='--')
        
    paths = [
        [
            (main_bracket_x + main_brackset_sign * (bracket_offset + bracket_width), inset_y_lim[1], main_ax),
            (main_bracket_x + main_brackset_sign * bracket_offset, inset_y_lim[1], main_ax),
            (main_bracket_x + main_brackset_sign * bracket_offset, inset_y_lim[0], main_ax),
            (main_bracket_x + main_brackset_sign * (bracket_offset + bracket_width), inset_y_lim[0], main_ax),
        ],
        [
            (inset_bracket_x + inset_bracket_sign * (bracket_offset + bracket_width), inset_y_lim[1], inset_ax),
            (inset_bracket_x + inset_bracket_sign * bracket_offset, inset_y_lim[1], inset_ax),
            (inset_bracket_x + inset_bracket_sign * bracket_offset, inset_y_lim[0], inset_ax),
            (inset_bracket_x + inset_bracket_sign * (bracket_offset + bracket_width), inset_y_lim[0], inset_ax),
        ],
    ]
    
    for path in paths:
        draw_line(path, linewidth=2)
    
    return inset_ax
