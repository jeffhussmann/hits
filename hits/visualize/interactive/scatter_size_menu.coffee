choice = cb_obj.value

if choice == ''
    choice = '_uniform_size'

scatter_source.data['_size'] = scatter_source.data[choice]

scatter_source.change.emit()