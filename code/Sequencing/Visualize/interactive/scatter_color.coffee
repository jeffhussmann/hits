models = cb_obj.document._all_models_by_name._dict

choice = cb_obj.value

if choice == ''
    models['scatter_source'].data['_color'] = models['scatter_source'].data['_black']
    models['scatter_source'].data['_selection_color'] = models['scatter_source'].data['_orange']
else
    models['scatter_source'].data['_color'] = models['scatter_source'].data[choice]
    models['scatter_source'].data['_selection_color'] = models['scatter_source'].data[choice]



models['scatter_source'].trigger('change')
