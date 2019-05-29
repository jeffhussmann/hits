models = cb_obj.document._all_models_by_name._dict

scatter_data = models['scatter_source'].data
label_data = models['filtered_source'].data
hist_data = models['histogram_source'].data

for axis in ['x', 'y']
    name = models[axis + '_0_menu'].value + ' ' + models[axis + '_1_menu'].value 

    scatter_data[axis] = scatter_data[name]
    label_data[axis] = label_data[name]

    for suffix in ['_all', '_bins_left', '_bins_right']
        hist_data[axis + suffix] = hist_data[name + suffix]

    models[axis + '_axis'].axis_label = name

# Call to recompute selection histograms.
models['scatter_selection_callback'].func(models['scatter_source'].selected, cb_data, require, exports)

models['scatter_source'].change.emit()
models['filtered_source'].change.emit()
models['histogram_source'].change.emit()