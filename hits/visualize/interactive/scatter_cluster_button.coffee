formatters =
    'clustered': {clustered_formatter}
    'original': {original_formatter}

cluster = not dendrogram_lines.visible
dendrogram_lines.visible = cluster

if cluster
    order_key = 'clustered'
else
    order_key = 'original'

data = heatmap_source.data

for k in ['r', 'x', 'x_name', 'y', 'y_name', 'color']
    data[k] = data[k + '_' + order_key]

for axis in [heatmap_x_axis, heatmap_y_axis]
    # Note: axis is a splattable list
    axis[0].formatter.code = formatters[order_key]

# Determine the new selection from the scatter axis labels.
# Note: each axis is a splattable list
x_name = x_axis[0].axis_label
y_name = y_axis[0].axis_label
num_pairs = data['x_name'].length
x_names = data['x_name']
y_names = data['y_name']
index = (i for i in [0..num_pairs] when x_names[i] == x_name and y_names[i] == y_name)[0]

heatmap_source.selected.indices = [index]

heatmap_source.change.emit()