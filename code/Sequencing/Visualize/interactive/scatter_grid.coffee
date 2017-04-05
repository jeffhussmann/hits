models = cb_obj.document._all_models_by_name._dict

choice = cb_obj.labels[cb_obj.active]

grids = (v for k, v of models when k == 'grid')[0]
for grid in grids
    grid.visible = choice == 'grid'

diagonals = (v for k, v of models when k == 'diagonal')[0]

for diagonal in diagonals
    diagonal.visible = choice == 'diagonal'
