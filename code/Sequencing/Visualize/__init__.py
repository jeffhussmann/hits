import igv_colors
import itertools
import scipy.stats
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def optional_ax(original_function):
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

def add_commas_to_yticks(ax):
    def commas_formatter(x, pos):
        return '{0:,}'.format(int(x))
    tick_formatter = matplotlib.ticker.FuncFormatter(commas_formatter)
    ax.yaxis.set_major_formatter(tick_formatter)

def enhanced_scatter(xs, ys, ax,
                     color_by_density=True,
                     do_fit=True,
                     show_p_value=True,
                     show_beta=True,
                     hists_height=0,
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
                     color_list=None,
                    ):
    xs = np.asarray(xs)
    ys = np.asarray(ys)

    same_lists = np.allclose(xs, ys)

    if color_by_density and not same_lists and len(xs) > 2:
        indices = np.arange(len(xs))
        np.random.shuffle(indices)
        random_indices = indices[:1000]

        sampled_points = np.vstack([xs[random_indices], ys[random_indices]])
        points = np.vstack([xs, ys])
        kernel = scipy.stats.gaussian_kde(sampled_points)
        colors = kernel(points)
    elif color_list:
        if isinstance(color_list, str):
            colors = color_list
        else:
            colors = np.asarray(color_list)
    else:
        colors = 'black'

    if same_lists:
        do_fit = False

    kwargs = {'cmap': matplotlib.cm.jet,
              's': marker_size,
              'linewidths' : (0.1,),
             }

    ax.scatter(xs, ys, c=colors, **kwargs)

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

    text_kwargs = {'xy': (x, y),
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

    if in_log_space:
        # Values have been log'ed before being passed in, so exponentiate them
        # to recover correct correlation.
        r, p = scipy.stats.pearsonr(2**xs, 2**ys)
    else:
        r, p = scipy.stats.pearsonr(xs, ys)

    if variance:
        r_part = '$r^2 = {:0.{digits}f}$'.format(r**2, digits=r_digits)
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
                **text_kwargs)

    original_xlims = ax.get_xlim()
    original_ylims = ax.get_ylim()

    if hists_height > 0:
        ax_x = ax.twinx()
        ax_x.hist(xs, alpha=0.9, histtype='step', bins=100)
        y_min, y_max = ax_x.get_ylim()
        ax_x.set_ylim(ymax=y_max / hists_height)

        ax_y = ax.twiny()
        ax_y.hist(ys, alpha=0.9, histtype='step', bins=100, orientation='horizontal')
        x_min, x_max = ax_y.get_xlim()
        ax_y.set_xlim(xmax=x_max / hists_height)

        ax_x.set_yticks([])
        ax_y.set_xticks([])

        ax.set_xlim(original_xlims)
        ax.set_ylim(original_ylims)
        
def draw_diagonal(ax, anti=False, color='black', **kwargs):
    if anti:
        xs, ys = [0, 1], [1, 0]
    else:
        xs, ys = [0, 1], [0, 1]
        if ax.get_xlim() != ax.get_ylim():
            raise ValueError('diagonal in non-equal axes')
    
    ax.plot(xs, ys,
            transform=ax.transAxes,
            color=color,
            **kwargs)

def label_scatter_plot(ax, xs, ys, labels, to_label,
                       vector='orthogonal',
                       initial_distance=50,
                       arrow_alpha=0.2,
                       manual_ratios=None,
                       manual_alignments=None,
                       text_kwargs={'size': 10},
                      ):
    def attempt_text(x, y, site, distance):
        if vector == 'orthogonal':
            x_offset = np.sign(x - y) * distance
            y_offset = -np.sign(x - y) * distance
            ha = 'center'
            va = 'top' if y_offset < 0 else 'bottom'
        elif vector == 'radial':
            norm = np.linalg.norm([x, y])
            x_offset = x * distance / norm
            y_offset = y * distance / norm
            ha, va = manual_alignments
        elif vector == 'sideways':
            x_offset = distance
            y_offset = 0
            ha, va = 'center', 'top'
        elif vector == 'manual':
            x_ratio, y_ratio = manual_ratios
            ha, va = manual_alignments
            x_offset = distance * x_ratio
            y_offset = distance * y_ratio

        text = ax.annotate(site,
                           xy=(x, y),
                           xycoords=('data', 'data'),
                           xytext=(x_offset, y_offset),
                           textcoords='offset points',
                           ha=ha,
                           va=va,
                           **text_kwargs)
        ax.figure.canvas.draw()

        return text, text.get_window_extent(), (x, y, x_offset, y_offset)

    ax.figure.canvas.draw()
    starting_labels = [ax.xaxis.get_label(), ax.yaxis.get_label()] + ax.get_yticklabels() + ax.get_xticklabels()
    bboxes = [label.get_window_extent() for label in starting_labels]

    tuples = itertools.izip(np.asarray(xs)[to_label],
                            np.asarray(ys)[to_label],
                            np.asarray(labels)[to_label],
                            )
    for x, y, label in tuples:
        distance = initial_distance
        text, bbox, coords = attempt_text(x, y, label, distance)
        while any(bbox.fully_overlaps(other_bbox) for other_bbox in bboxes):
            text.remove()
            distance += 10
            text, bbox, coords = attempt_text(x, y, label, distance)
            if distance >= 500:
                break
        
        x, y, x_offset, y_offset = coords
        ax.annotate('',
                    xy=(x, y),
                    xycoords=('data', 'data'),
                    xytext=(x_offset, y_offset),
                    textcoords=('offset points', 'offset points'),
                    arrowprops={'arrowstyle': '->', 'alpha': 0.5},
                   )

        bboxes.append(bbox)

def evenly_spaced_jet_colors(n):
    jet_colors = [matplotlib.cm.jet(x) for x in np.linspace(0, 1, n)]
    return jet_colors

def assign_colors_by_ranks(key_to_ranks, by_last=False):
    colors = evenly_spaced_jet_colors(len(key_to_ranks))
    if by_last:
        index = -1
    else:
        index = 0
    key_to_color = {key: colors[key_to_ranks[key][index]] for key in key_to_ranks}
    return key_to_color
