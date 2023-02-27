column_names = {column_names}

is_case_sensitive = case_sensitive.active.length > 0
if not is_case_sensitive
    queries = search_input.value.toLowerCase().split('|')
    possibly_lowercase = (t) -> t.toString().toLowerCase()
else
    queries = search_input.value.split('|')
    possibly_lowercase = (t) -> t

matches_any_query = (t) -> 
    matches_query = (query) -> possibly_lowercase(t).indexOf(query) > -1 
    queries.some(matches_query)

all_matches = []
if search_input.value != ''
    for column in column_names
        targets = scatter_source.data[column]
        matches_from_column = (i for t, i in targets when matches_any_query(t) and i not in all_matches)
        all_matches.push matches_from_column...

scatter_source.selected.indices = all_matches