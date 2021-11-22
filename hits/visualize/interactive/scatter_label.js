var choice;

choice = cb_obj.value;

scatter_source.data['_label'] = scatter_source.data[choice];

filtered_source.data['_label'] = filtered_source.data[choice];

filtered_source.change.emit();
