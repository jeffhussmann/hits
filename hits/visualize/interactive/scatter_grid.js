var axes_line, axes_lines, choice, diagonal, diagonals, grid, i, j, k, len, len1, len2;

choice = cb_obj.labels[cb_obj.active];

for (i = 0, len = fig_grid.length; i < len; i++) {
  grid = fig_grid[i];
  grid.visible = choice === 'grid';
}

diagonals = fig_diagonal;

if (!Array.isArray(diagonals)) {
  diagonals = [diagonals];
}

for (j = 0, len1 = diagonals.length; j < len1; j++) {
  diagonal = diagonals[j];
  diagonal.visible = choice.includes('diagonal');
}

axes_lines = fig_axes_line;

if (!Array.isArray(axes_lines)) {
  axes_lines = [axes_lines];
}

for (k = 0, len2 = axes_lines.length; k < len2; k++) {
  axes_line = axes_lines[k];
  axes_line.visible = choice.includes('axes');
}
