log_scale = {log_scale}
identical_bins = {identical_bins}

ranges =
    'x': x_range
    'y': y_range
    'hist_x': hist_x_range
    'hist_y': hist_y_range 

xs = scatter_source.data.x
ys = scatter_source.data.y

clean_xs = []
clean_ys = []

if log_scale
    number_is_clean = (x) -> x isnt 0 and not Number.isNaN(x) and x isnt 'NaN'
else
    number_is_clean = (x) -> not Number.isNaN(x) and x isnt 'NaN'

index_is_clean = (i) -> number_is_clean(xs[i]) and number_is_clean(ys[i])

for i in [0...xs.length]
    if index_is_clean(i)
        clean_xs.push(xs[i])
        clean_ys.push(ys[i])

clean_values =
    'x': clean_xs
    'y': clean_ys

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

    # If identical_bins, want ranges of each histogram to be
    # the same so that heights can be meaningfully compared.
    # If not, maximize dynamic range. 
    if not identical_bins
        hist_counts = histogram_source.data[axis_name + '_all']
        hist_max = Math.max(hist_counts...)
        ranges['hist_' + axis_name].end = hist_max
