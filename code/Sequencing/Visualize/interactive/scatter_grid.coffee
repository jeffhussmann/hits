models = cb_obj.document._all_models_by_name._dict

grids = (v for k, v of models when k == 'grid')[0]
for grid in grids
    grid.visible = not grid.visible

diagonals = (v for k, v of models when k == 'diagonal')[0]
for diagonal in diagonals
    diagonal.visible = not diagonal.visible
