num_selected = cb_obj.indices.length

axes = 
    'x': x_axis
    'y': y_axis

if num_selected == 0
    # Selection was cleared with ESC. Recover the old selection from the axis
    # labels.
    x_name = axes['x'].axis_label
    y_name = axes['y'].axis_label
    num_pairs = heatmap_source.data['x_name'].length
    x_names = heatmap_source.data['x_name']
    y_names = heatmap_source.data['y_name']
    index = (i for i in [0..num_pairs] when x_names[i] == x_name and y_names[i] == y_name)[0]
else
    # In case of multiple selection with shift key, only keep the most recent.
    index = cb_obj.indices[num_selected - 1]

cb_obj.indices = [index]

for axis in ['x', 'y']
    name = heatmap_source.data[axis + '_name'][index]

    scatter_source.data[axis] = scatter_source.data[name]
    filtered_source.data[axis] = filtered_source.data[name]

    for suffix in ['_all', '_bins_left', '_bins_right']
        histogram_source.data[axis + suffix] = histogram_source.data[name + suffix]
    
    axes[axis].axis_label = name

scatter_source.change.emit()
filtered_source.change.emit()
histogram_source.change.emit()