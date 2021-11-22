squeeze = (possibly_array) ->
    if Array.isArray(possibly_array)
        squeezed = possibly_array[0]
    else
        squeezed = possibly_array
    return squeezed
    
choice = squeeze cb_obj.value

if choice == ''
    main_key = '_black'
    selection_key = '_orange'
else
    main_key = choice
    selection_key = choice

scatter_source.data['_color'] = scatter_source.data[main_key]
scatter_source.data['_selection_color'] = scatter_source.data[selection_key]

scatter_source.change.emit()