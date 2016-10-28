scatter_data = scatter_source.data
label_data = label_source.data

x_name = x_menu.value
y_name = y_menu.value

scatter_data.x = scatter_data[x_name]
scatter_data.y = scatter_data[y_name]

label_data.x = label_data[x_name]
label_data.y = label_data[y_name]

xaxis.axis_label = x_name
yaxis.axis_label = y_name

scatter_source.trigger('change')
label_source.trigger('change')
