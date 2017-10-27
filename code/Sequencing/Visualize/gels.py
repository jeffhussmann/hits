import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image
import os
import scipy.signal
import numpy as np
import ipywidgets
from collections import defaultdict

def load_labels(image_fn):
    labels_fn = '{0}_labels.txt'.format(image_fn[:-len('.tif')])
    labels = [line.strip() for line in open(labels_fn)]
    return labels

def identify_lane_boundaries(image):
    cols = image.sum(axis=0)
    peaks, = scipy.signal.argrelmax(sum_over_window(cols, 20))
    mean_width = np.mean(np.diff(peaks))
    either_side = int(mean_width * 0.4)
    boundaries = [(p - either_side, p + either_side) for p in peaks]
    return boundaries

def extract_profiles(image, boundaries):
    profiles = []
    for start, end in boundaries:
        lane = image[:, start:end + 1]
        profile = lane.sum(axis=1)[::-1]
        profiles.append(profile)
    return profiles

def sum_over_window(array, either_side):
    summed = np.zeros(len(array))
    cumulative = np.concatenate(([0], np.cumsum(array)))
    window = 2 * either_side + 1
    summed[either_side:-either_side] = cumulative[window:] - cumulative[:-window]
    return summed

def get_ladder100_peaks(ys):
    ys = np.asarray(ys)
    lengths = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1517]
    
    relative_maxes,  = scipy.signal.argrelmax(ys, order=3)
    threshold = np.percentile(ys, 75)
    candidates = relative_maxes[ys[relative_maxes] > threshold]
    
    xs = np.arange(len(ys))[candidates]
    if len(xs) == len(lengths):
        peaks = dict(zip(lengths, xs))
    else:
        # Peaks are probably too close together to call.
        # Give up on the higher ones.
        peaks = dict(zip(lengths[:5], xs[:5]))
        peaks[lengths[-1]] = xs[-1]
        
    return peaks

def get_ladder10_peaks(ys, ladder100_peaks):
    # Strategy: identify the 100 bp beak by finding the relative max closest
    # to the 100 bp peak in the 100bp ladder.
    # March up and down relative maxes from there to assign 10 bp peaks.
    # To avoid false positive peaks between the lower ladder components,
    # require them to be the max in a larger index winow (via the order
    # argument to argrelmax).
    
    ys = np.asarray(ys)
    
    def get_xs(order):
        relative_maxs, = scipy.signal.argrelmax(ys, order=order)
        xs = np.arange(len(ys))[relative_maxs]
        distance = np.abs(xs - ladder100_peaks[100])
        index_100 = np.argmin(distance)
        return xs, index_100
        
    # below 100
    xs, index_100 = get_xs(10)
    peaks = {100: xs[index_100]}
    
    for length, i in zip(np.arange(90, 10, -10), range(index_100 - 1, 0, -1)):
        peaks[length] = xs[i]
    
    # above 100
    xs, index_100 = get_xs(3)
    
    for length, i in zip(np.arange(110, 160, 10), range(index_100 + 1, len(xs))):
        peaks[length] = xs[i]
        
    return peaks

def plot_gel(image_fn, vertical_range=(0, 1), show_ladders=False, highlight=None):
    image = matplotlib.image.imread(image_fn)
    boundaries = identify_lane_boundaries(image)
    profiles = extract_profiles(image, boundaries)
    
    labels = load_labels(image_fn)
    
    lanes = dict(zip(labels, profiles))
    
    if len(labels) != len(boundaries):
        raise ValueError('mismatch between number of labels and number of detected lanes')
        
    rows, cols = image.shape
    im_width = 8. * cols / rows

    gridspec_kw = dict(width_ratios=[16, im_width], wspace=0)
    
    fig, (line_ax, im_ax) = plt.subplots(1, 2, figsize=(16 + im_width, 8), gridspec_kw=gridspec_kw)
    
    kwargs = {
        'line': {
            'highlight': dict(alpha=0.95, linewidth=1),
            'nonhighlight': dict(alpha=0.1, linewidth=0.5),
            'uniform': dict(alpha=0.9, linewidth=0.7),
        },
        'boundary': {
            'highlight': dict(alpha=0.95, linewidth=1.5),
            'nonhighlight': dict(alpha=0.2, linewidth=1),
            'uniform': dict(alpha=0.9, linewidth=1),
        },
        'text': {
            'highlight': dict(alpha=0.95),
            'nonhighlight': dict(alpha=0.1),
            'uniform': dict(alpha=0.9),
        },
    }
    
    label_to_kwargs = defaultdict(dict)

    if highlight is None:
        highlight = []
    
    for i, label in enumerate(labels):
        color = 'C{0}'.format(i)
        
        if len(highlight) > 0:
            if label in highlight:
                key = 'highlight'
            else:
                key = 'nonhighlight'
        else:
            key = 'uniform'
        
        for kind in kwargs:
            copy = kwargs[kind][key].copy()
            copy['color'] = color
            label_to_kwargs[kind][label] = copy
    
    for i, label in enumerate(labels):
        if 'ladder' in label and not show_ladders:
            continue
            
        ys = lanes[label]
        xs = np.arange(len(ys))
        
        line_ax.plot(xs, ys, 'o-', label=label, markersize=0.3, **label_to_kwargs['line'][label])

    ladder100_peaks = get_ladder100_peaks(lanes['ladder100'])

    ladder10_peaks = get_ladder10_peaks(lanes['ladder10'], ladder100_peaks)
    
    # Only include the 100 bp peak from the 10 bp ladder.
    ladder100_peaks.pop(100)

    peaks = ladder100_peaks.items() + ladder10_peaks.items()

    major_peaks = [
        100,
        200,
        500,
        1000,
        1517,
    ]

    for length, x in peaks:
        alpha = 0.3 if length in major_peaks else 0.05
        line_ax.axvline(x, color='black', alpha=alpha)

    line_ax.set_yticks([])
    
    line_ax.set_xticks([x for length, x in peaks])
    line_ax.set_xticklabels([str(length) for length, x in peaks], rotation=-90, ha='center', size=8)
    
    legend = line_ax.legend(loc='upper left', framealpha=0.5)
    
    im_ax.imshow(image, cmap=matplotlib.cm.binary_r)
    im_ax.set_xticks([])
    im_ax.set_yticks([])
    
    for i, (label, (start, end)) in enumerate(zip(labels, boundaries)):
        im_ax.axvline(start, **label_to_kwargs['boundary'][label])
        im_ax.axvline(end, **label_to_kwargs['boundary'][label])
        im_ax.annotate(label,
                       xy=(np.mean((start, end)), 1),
                       xytext=(0, 2),
                       xycoords=('data', 'axes fraction'),
                       textcoords='offset points',
                       rotation=45,
                       va='bottom',
                       **label_to_kwargs['text'][label])
        
    for text in legend.get_texts():
        label = text.get_text()
        text.set(**label_to_kwargs['text'][label])
   
    x_min, x_max = map(int, np.asarray(vertical_range) * len(xs))
    line_ax.set_xlim(x_min, x_max)
    
    im_ax.autoscale(False)

    if vertical_range != (0, 1):
        y_min, y_max = vertical_range

        if y_min == 0:
            y_min = 0.002
        if y_max == 1:
            y_max= 0.999

        line_kwargs = dict(transform=im_ax.transAxes, color='white')
        im_ax.plot([0.005, 0.005, 0.05], [y_min + 0.05, y_min, y_min], **line_kwargs)
        im_ax.plot([0.005, 0.005, 0.05], [y_max - 0.05, y_max, y_max], **line_kwargs)
        im_ax.plot([1 - 0.005, 1 - 0.005, 1 - 0.05], [y_min + 0.05, y_min, y_min], **line_kwargs)
        im_ax.plot([1 - 0.005, 1 - 0.005, 1 - 0.05], [y_max - 0.05, y_max, y_max], **line_kwargs)
        
    head, tail = os.path.split(image_fn)
    line_ax.set_title(tail)
        
    return fig

def plot_gels_interactive(image_fns, **kwargs):

    def make_tab(image_fn):
        def generate_figure(highlight, vertical_range):
            fig = plot_gel(image_fn, highlight=highlight, vertical_range=vertical_range)
            plt.show()
            return fig

        widgets = {
            'vertical_range': ipywidgets.FloatRangeSlider(value=[0, 1],
                                                          continuous_update=False,
                                                          min=0,
                                                          max=1,
                                                          step=0.01,
                                                          layout=ipywidgets.Layout(height='200px'),
                                                          style={'description_width': 'initial'},
                                                          orientation='vertical',
                                                         ),
            'highlight': ipywidgets.SelectMultiple(options=load_labels(image_fn),
                                                   value=[],
                                                   layout=ipywidgets.Layout(height='200px'),
                                                  ),
            'save': ipywidgets.Button(description='Save'),
            'file_name': ipywidgets.Text(value=os.environ['HOME'] + '/', description='file name'),
        }

        def save(button):
            fig = interactive.result
            fn = widgets['file_name'].value
            fig.savefig(fn, bbox_inches='tight')

        widgets['save'].on_click(save)

        interactive = ipywidgets.interactive(generate_figure, **widgets)
        output = interactive.children[-1]
        output.layout.height = '870px'
        interactive.update()

        rows = [
            [output],
            [widgets['highlight'], widgets['vertical_range'], widgets['file_name'], widgets['save']],
        ]
        cols = ipywidgets.VBox([ipywidgets.HBox(row) for row in rows])
        return cols 

    tab = ipywidgets.Tab()
    tab.children = [make_tab(image_fn) for image_fn in image_fns]
    for i, image_fn in enumerate(image_fns):
        head, tail = os.path.split(image_fn)
        tab.set_title(i, tail)

    return tab
