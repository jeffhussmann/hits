for axis_name in ['x', 'y']
    if axis_name == 'x'
        menu_0 = x_menu_0
        menu_1 = x_menu_1
        axis_title = x_axis_title
        axis_subtitle = x_axis_subtitle
    else
        menu_0 = y_menu_0
        menu_1 = y_menu_1
        axis_title = y_axis_title
        axis_subtitle = y_axis_subtitle

    title = menu_0.value[0]
    subtitle = menu_1.value[0]
    name = title + ' ' + subtitle 

    scatter_source.data[axis_name] = scatter_source.data[name]
    filtered_source.data[axis_name] = filtered_source.data[name]

    for suffix in ['_all', '_bins_left', '_bins_right']
        histogram_source.data[axis_name + suffix] = histogram_source.data[name + suffix]

    axis_title.text = title
    axis_subtitle.text = subtitle

# Note: need to re-add calculation of histograms.

scatter_source.change.emit()
filtered_source.change.emit()
histogram_source.change.emit()