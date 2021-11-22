column_names = {column_names}

is_case_sensitive = case_sensitive.active.length > 0
if not is_case_sensitive
    query = search_input.value.toLowerCase()
    possibly_lowercase = (t) -> t.toString().toLowerCase()
else
    query = cb_obj.value
    possibly_lowercase = (t) -> t

all_matches = []
if query != ''
    for column in column_names
        targets = scatter_source.data[column]
        matches = (i for t, i in targets when possibly_lowercase(t).indexOf(query) > -1 and i not in all_matches)
        all_matches.push matches...

scatter_source.selected.indices = all_matches