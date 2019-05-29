models = cb_obj.document._all_models_by_name._dict

scatter_data = models['scatter_source'].data
label_data = models['filtered_source'].data
hist_data = models['histogram_source'].data
heatmap_data = models['heatmap_source'].data

num_selected = cb_obj.indices.length

if num_selected == 0
    # Selection was cleared with ESC. Recover the old selection from the axis
    # labels.
    x_name = models['x_axis'].axis_label
    y_name = models['y_axis'].axis_label
    num_pairs = heatmap_data['x_name'].length
    x_names = heatmap_data['x_name']
    y_names = heatmap_data['y_name']
    index = (i for i in [0..num_pairs] when x_names[i] == x_name and y_names[i] == y_name)[0]
else
    # In case of multiple selection with shift key, only keep the most recent.
    index = cb_obj.indices[num_selected - 1]

cb_obj.indices = [index]

for axis in ['x', 'y']
    name = heatmap_data[axis + '_name'][index]

    scatter_data[axis] = scatter_data[name]
    label_data[axis] = label_data[name]

    for suffix in ['_all', '_bins_left', '_bins_right']
        hist_data[axis + suffix] = hist_data[name + suffix]
    
    models[axis + '_axis'].axis_label = name

# Call to recompute selection histograms.
models['scatter_selection_callback'].func(models['scatter_source'].selected, 'from_heatmap', require, exports)

models['scatter_source'].change.emit()
models['filtered_source'].change.emit()
models['histogram_source'].change.emit()