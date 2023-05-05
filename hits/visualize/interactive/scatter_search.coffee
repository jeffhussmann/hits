column_names = {column_names}

search_input = search_input.value

is_case_sensitive = case_sensitive.active.length > 0
if not is_case_sensitive
    search_input = search_input.toLowerCase()
    possibly_lowercase = (t) -> t.toString().toLowerCase()
else
    possibly_lowercase = (t) -> t

matches_query = (target) ->
    (query) -> possibly_lowercase(target).indexOf(query) > -1 

if search_input.indexOf('|') > -1
    queries = search_input.split('|')

    matches_relevant_query = (target) -> queries.some(matches_query(target))

else if search_input.indexOf('&') > -1
    queries = search_input.split('&')

    matches_relevant_query = (target) -> queries.every(matches_query(target))

else
    matches_relevant_query = (target) -> matches_query(target)(search_input)


all_matches = []
if search_input.value != ''
    for column in column_names
        targets = scatter_source.data[column]
        matches_from_column = (i for target, i in targets when matches_relevant_query(target) and i not in all_matches)
        all_matches.push matches_from_column...

scatter_source.selected.indices = all_matches