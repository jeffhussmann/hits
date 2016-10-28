full_data = source.data
filtered_data = table.source.data
indices = cb_obj.selected['1d'].indices

for key, values of full_data
    filtered_data[key] = (values[i] for i in indices)

table.trigger('change')
labels.trigger('change')
