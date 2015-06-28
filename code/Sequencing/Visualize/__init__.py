import igv_colors
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
                    ):
    same_lists = np.allclose(xs, ys)

    if color_by_density and not same_lists and len(xs) > 2:
        indices = np.arange(len(xs))
        np.random.shuffle(indices)
        random_indices = indices[:1000]

        sampled_points = np.vstack([xs[random_indices], ys[random_indices]])
        points = np.vstack([xs, ys])
        kernel = scipy.stats.gaussian_kde(sampled_points)
        colors = kernel(points)
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
        ax.plot(x_lims, fit_function(x_lims), **fit_line_kwargs)
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

    if hists_height > 0:
        ax_x = ax.twinx()
        ax_x.hist(xs, alpha=0.3, histtype='step', bins=100)
        y_min, y_max = ax_x.get_ylim()
        ax_x.set_ylim(ymax=y_max / hists_height)

        ax_y = ax.twiny()
        ax_y.hist(ys, alpha=0.3, histtype='step', bins=100, orientation='horizontal')
        x_min, x_max = ax_y.get_xlim()
        ax_y.set_xlim(xmax=x_max / hists_height)

        ax_x.set_yticks([])
        ax_y.set_xticks([])
        
def draw_diagonal(ax, color='black', **kwargs):
    if ax.get_xlim() != ax.get_ylim():
        raise ValueError('diagonal in non-equal axes')

    ax.plot([0, 1], [0, 1],
            transform=ax.transAxes,
            color=color,
            **kwargs)
