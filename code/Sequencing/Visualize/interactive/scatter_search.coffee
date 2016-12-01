models = cb_obj.document._all_models_by_name._dict

query = cb_obj.value

all_matches = []
if query != ''
    for column in {columns}
        targets = models['scatter_source'].data[column]
        matches = (i for t, i in targets when t.indexOf(query) > -1 and i not in all_matches)
        all_matches.push matches...

models['scatter_source'].selected['1d'].indices = all_matches
models['scatter_source'].callback.func(models['scatter_source'])
