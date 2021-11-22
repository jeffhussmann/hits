choice = cb_obj.labels[cb_obj.active]

for grid in fig_grid
    grid.visible = choice == 'grid'

diagonals = fig_diagonal
if not Array.isArray(diagonals)
    diagonals = [diagonals]

for diagonal in diagonals
    diagonal.visible = choice.includes('diagonal')

axes_lines = fig_axes_line
if not Array.isArray(axes_lines)
    axes_lines = [axes_lines]

for axes_line in axes_lines
    axes_line.visible = choice.includes('axes')