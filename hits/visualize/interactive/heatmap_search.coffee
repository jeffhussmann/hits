data = quad_source.data

row_labels = {row_labels}
col_labels = {col_labels}

query = cb_obj.value

if query != ''
    row_matches = (i for r, i in row_labels when r.indexOf(query) > -1)
    row_alphas = (0.3 for r in row_labels)
    row_alphas[i] = 1 for i in row_matches

    col_matches = (i for c, i in col_labels when c.indexOf(query) > -1)
    col_alphas = (0.3 for c in col_labels)
    col_alphas[i] = 1 for i in col_matches
else
    row_matches = []
    row_alphas = (1 for r in row_labels)

    col_matches = []
    col_alphas = (1 for c in col_labels)

data['left'] = [({lower_bound} for i in row_matches)..., ({lower_bound} + i for i in col_matches)...]
data['right'] = [({x_upper_bound} for i in row_matches)..., ({lower_bound} + i + 1 for i in col_matches)...]
data['bottom'] = [({y_upper_bound} - i - 1 for i in row_matches)..., ({lower_bound} for i in col_matches)...]
data['top'] = [({y_upper_bound} - i for i in row_matches)..., ({y_upper_bound} for i in col_matches)...]

row_label_source.data['alpha'] = row_alphas
col_label_source.data['alpha'] = col_alphas

quad_source.change.emit()
row_label_source.change.emit()
col_label_source.change.emit()