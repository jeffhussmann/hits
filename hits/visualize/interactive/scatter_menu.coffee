squeeze = (possibly_array) ->
    if Array.isArray(possibly_array)
        squeezed = possibly_array[0]
    else
        squeezed = possibly_array
    return squeezed

menus = 
    'x': x_menu
    'y': y_menu

axes = 
    'x': x_axis
    'y': y_axis

for axis_name in ['x', 'y']
    menu = menus[axis_name]
    axis = axes[axis_name]

    name = squeeze menu.value

    scatter_source.data[axis_name] = scatter_source.data[name]
    filtered_source.data[axis_name] = filtered_source.data[name]

    for suffix in ['_all', '_bins_left', '_bins_right']
        histogram_source.data[axis_name + suffix] = histogram_source.data[name + suffix]

    axis[0].axis_label = name

# Note: need to re-add calculation of histograms.

scatter_source.change.emit()
filtered_source.change.emit()
histogram_source.change.emit()