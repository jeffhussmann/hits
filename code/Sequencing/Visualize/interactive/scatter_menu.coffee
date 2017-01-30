models = cb_obj.document._all_models_by_name._dict

scatter_data = models['scatter_source'].data
label_data = models['labels_source'].data

x_name = models['x_menu'].value
y_name = models['y_menu'].value

scatter_data.x = scatter_data[x_name]
scatter_data.y = scatter_data[y_name]

label_data.x = label_data[x_name]
label_data.y = label_data[y_name]

xaxis.axis_label = x_name
yaxis.axis_label = y_name

models['scatter_source'].trigger('change')
models['labels_source'].trigger('change')
