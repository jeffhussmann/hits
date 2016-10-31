query = cb_obj.value

all_matches = []
if query != ''
    for column in {columns}
        targets = scatter_source.data[column]
        matches = (i for t, i in targets when t.indexOf(query) > -1 and i not in all_matches)
        all_matches.push matches...

scatter_source.selected['1d'].indices = all_matches

scatter_source.callback.func(labels, scatter_source, table, scatter_source)
