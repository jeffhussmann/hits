models = cb_obj.document._all_models_by_name._dict

query = cb_obj.value

selection = []
if query != ''
    bools = models['scatter_source'].data[query]
    selection = (i for b, i in bools when b)

models['scatter_source'].selected['1d'].indices = selection
models['scatter_source'].callback.func(models['scatter_source'])
