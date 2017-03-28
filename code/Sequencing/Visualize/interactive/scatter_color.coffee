models = cb_obj.document._all_models_by_name._dict

choice = cb_obj.value

if choice == ''
    choice = '_no_color'

models['scatter_source'].data['_color'] = models['scatter_source'].data[choice]
models['scatter_source'].trigger('change')
