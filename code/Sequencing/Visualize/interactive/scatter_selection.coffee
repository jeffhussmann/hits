models = cb_obj.document._all_models_by_name._dict

full_data = models['scatter_source'].data
filtered_data = models['labels_source'].data
indices = cb_obj.selected['1d'].indices

for key, values of full_data
    filtered_data[key] = (values[i] for i in indices)

if (models['table']?)
    models['table'].trigger('change')

models['labels_source'].trigger('change')
