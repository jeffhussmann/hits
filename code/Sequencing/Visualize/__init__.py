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
                     hists_height=0,
                     marker_size=4,
                     text_size=14,
                     text_location='lower right',
                     fit_line_kwargs={'color': 'black',
                                      'alpha': 0.5,
                                     },
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
        colors = np.ones_like(xs)

    if same_lists:
        do_fit = False

    kwargs = {'cmap': matplotlib.cm.jet,
              's': marker_size,
              'linewidths' : (0.1,),
             }

    ax.scatter(xs, ys, c=colors, **kwargs)

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
        y_offset = 15

    text_kwargs = {'xy': (x, y),
                   'xycoords': 'axes fraction',
                   'textcoords': 'offset points',
                   'horizontalalignment': horizontal_alignment,
                   'verticalalignment': vertical_alignment,
                   'fontsize': text_size,
                   'family': 'serif',
                  }
    
    if do_fit:
        fit = np.polyfit(xs, ys, 1)
        beta, _ = fit
        fit_function = np.poly1d(fit)
        x_lims = ax.get_xlim()
        ax.plot(x_lims, fit_function(x_lims), **fit_line_kwargs)
        ax.set_xlim(*x_lims)
        
        ax.annotate(r'$\beta$ = {:0.2f}'.format(beta),
                    xytext=(x_sign * x_offset, y_sign * y_offset * 2),
                    **text_kwargs)
    
    r, p = scipy.stats.pearsonr(xs, ys)
    if show_p_value:
        text = 'r = {:0.2f}, p={:0.2e}'.format(r, p)
    else:
        text = 'r = {:0.2f}'.format(r)

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
        
def draw_diagonal(ax, color='black', alpha=0.9, **kwargs):
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color=color, alpha=alpha, **kwargs)
