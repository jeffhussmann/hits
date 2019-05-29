models = cb_obj.document._all_models_by_name._dict

choice = cb_obj.labels[cb_obj.active]

for grid in models['grid']
    grid.visible = choice == 'grid'

diagonals = models['diagonal']
if not Array.isArray(diagonals)
    diagonals = [diagonals]

for diagonal in diagonals
    diagonal.visible = choice.includes('diagonal')

axes_lines = models['axes_line']
if not Array.isArray(axes_lines)
    axes_lines = [axes_lines]

for axes_line in axes_lines
    axes_line.visible = choice.includes('axes')