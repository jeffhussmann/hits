models = cb_obj.document._all_models_by_name._dict

lines = (v for k, v of models when k.startsWith('line_'))
for line in lines
    line.data_source.data['y'] = line.data_source.data[cb_obj.value]
    line.data_source.trigger('change')
