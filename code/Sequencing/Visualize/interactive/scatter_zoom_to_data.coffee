models = cb_obj.document._all_models_by_name._dict

log_scale = {log_scale}
identical_bins = {identical_bins}

scatter_data = models['scatter_source'].data
ranges = {{
    'x': models['x_range'],
    'y': models['y_range'],
}}

xs = scatter_data.x
ys = scatter_data.y

clean_xs = []
clean_ys = []

if log_scale
    is_clean = (i) -> xs[i] isnt 0 and xs[i] isnt 'NaN' and ys[i] isnt 0 and ys[i] isnt 'NaN'
else
    is_clean = (i) -> xs[i] isnt 'NaN' and ys[i] isnt 'NaN'

for i in [0...xs.length]
    if is_clean(i)
        clean_xs.push(xs[i])
        clean_ys.push(ys[i])

clean_values = {{
    'x': clean_xs,
    'y': clean_ys,
}}

for axis_name in ['x', 'y']
    min = Math.min(clean_values[axis_name]...)
    max = Math.max(clean_values[axis_name]...)
    
    if log_scale
        ranges[axis_name].start = min / 2
        ranges[axis_name].end = max * 2
    else
        buffer = (max - min) * 0.05
        ranges[axis_name].start = min - buffer
        ranges[axis_name].end = max + buffer

    if not identical_bins
        hist_counts = models['histogram_source'].data[axis_name + '_all']
        hist_max = Math.max(hist_counts...)
        hist_range = models['hist_' + axis_name + '_range']
        hist_range.end = hist_max
