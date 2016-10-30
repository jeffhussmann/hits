models = cb_obj.document._all_models_by_name._dict

sources = (v for k, v of models when k.startsWith('source_'))

for source in sources
    source.data['y'] = source.data[cb_obj.value]
    source.trigger('change')
