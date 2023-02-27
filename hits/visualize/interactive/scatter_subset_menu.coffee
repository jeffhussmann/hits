subset_indices = {subset_indices}

query = cb_obj.value

selection = []
if query != ''
    selection = subset_indices[query]

scatter_source.selected.indices = selection