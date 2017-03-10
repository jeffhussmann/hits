models = cb_obj.document._all_models_by_name._dict

scatter_data = models['scatter_source'].data
label_data = models['labels_source'].data
hist_data = models['histogram_source'].data

squeeze = (possibly_array) ->
    if Array.isArray(possibly_array)
        squeezed = possibly_array[0]
    else
        squeezed = possibly_array
    return squeezed

x_name = squeeze models['x_menu'].value
y_name = squeeze models['y_menu'].value

scatter_data.x = scatter_data[x_name]
scatter_data.y = scatter_data[y_name]

label_data.x = label_data[x_name]
label_data.y = label_data[y_name]

hist_data['x_all'] = hist_data[x_name]
hist_data['y_all'] = hist_data[y_name]

models['x_axis'].axis_label = x_name
models['y_axis'].axis_label = y_name

# Call to recompute selection histograms.
models['scatter_source'].callback.func(models['scatter_source'])

models['scatter_source'].trigger('change')
models['labels_source'].trigger('change')
models['histogram_source'].trigger('change')
