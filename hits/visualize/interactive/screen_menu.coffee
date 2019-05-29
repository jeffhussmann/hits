models = cb_obj.document._all_models_by_name._dict

scatter_data = models['scatter_source'].data
label_data = models['filtered_source'].data

dataset_name = models['dataset_menu'].value[0]
outcome_name = models['outcome_menu'].value[0]

nt_fraction_values = {nt_fractions}

for key in ['frequency', 'ys', 'gene_p_up', 'gene_p_down', 'total_UMIs']
    full_name = dataset_name + '_' + outcome_name + '_' + key
    scatter_data[key] = scatter_data[full_name]
    label_data[key] = label_data[full_name]

full_name = dataset_name + '_' + outcome_name
models['nt_fraction'].location = nt_fraction_values[full_name]

models['title'].text = dataset_name + '     ' + outcome_name

models['title'].change.emit()
models['nt_fraction'].change.emit()
models['scatter_source'].change.emit()
models['filtered_source'].change.emit()
