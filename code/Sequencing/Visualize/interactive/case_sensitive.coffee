models = cb_obj.document._all_models_by_name._dict
models['search'].callback.func(models['search'])
