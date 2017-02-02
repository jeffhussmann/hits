models = cb_obj.document._all_models_by_name._dict
console.log models

scatter_data = models['scatter_source'].data
label_data = models['labels_source'].data

num_selected = cb_obj.selected['1d'].indices.length

if num_selected == 0
    # Selection was cleared with ESC. Recover the old selection from the axis
    # labels.
    x_name = models['x_axis'].axis_label
    y_name = models['y_axis'].axis_label
    num_pairs = cb_obj.data['x_name'].length
    x_names = cb_obj.data['x_name']
    y_names = cb_obj.data['y_name']
    index = (i for i in [0..num_pairs] when x_names[i] == x_name and y_names[i] == y_name)[0]
else
    # In case of multiple selection with shift key, only keep the most recent.
    index = cb_obj.selected['1d'].indices[num_selected - 1]

cb_obj.selected['1d'].indices = [index]
x_name = cb_obj.data['x_name'][index]
y_name = cb_obj.data['y_name'][index]

scatter_data.x = scatter_data[x_name]
scatter_data.y = scatter_data[y_name]

label_data.x = label_data[x_name]
label_data.y = label_data[y_name]

models['x_axis'].axis_label = x_name
models['y_axis'].axis_label = y_name

models['scatter_source'].trigger('change')
models['labels_source'].trigger('change')
