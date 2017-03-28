models = cb_obj.document._all_models_by_name._dict

choice = cb_obj.value

models['scatter_source'].data['_label'] = models['scatter_source'].data[choice]
models['labels_source'].data['_label'] = models['labels_source'].data[choice]
models['labels_source'].trigger('change')
