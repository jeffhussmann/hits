var axes, axis, i, index, j, k, len, len1, name, num_pairs, num_selected, ref, ref1, suffix, x_name, x_names, y_name, y_names;

num_selected = cb_obj.indices.length;

axes = {
  'x': x_axis,
  'y': y_axis
};

if (num_selected === 0) {
  // Selection was cleared with ESC. Recover the old selection from the axis
  // labels.
  x_name = axes['x'].axis_label;
  y_name = axes['y'].axis_label;
  num_pairs = heatmap_source.data['x_name'].length;
  x_names = heatmap_source.data['x_name'];
  y_names = heatmap_source.data['y_name'];
  index = ((function() {
    var j, ref, results;
    results = [];
    for (i = j = 0, ref = num_pairs; (0 <= ref ? j <= ref : j >= ref); i = 0 <= ref ? ++j : --j) {
      if (x_names[i] === x_name && y_names[i] === y_name) {
        results.push(i);
      }
    }
    return results;
  })())[0];
} else {
  // In case of multiple selection with shift key, only keep the most recent.
  index = cb_obj.indices[num_selected - 1];
}

cb_obj.indices = [index];

ref = ['x', 'y'];
for (j = 0, len = ref.length; j < len; j++) {
  axis = ref[j];
  name = heatmap_source.data[axis + '_name'][index];
  scatter_source.data[axis] = scatter_source.data[name];
  filtered_source.data[axis] = filtered_source.data[name];
  ref1 = ['_all', '_bins_left', '_bins_right'];
  for (k = 0, len1 = ref1.length; k < len1; k++) {
    suffix = ref1[k];
    histogram_source.data[axis + suffix] = histogram_source.data[name + suffix];
  }
  axes[axis].axis_label = name;
}

scatter_source.change.emit();

filtered_source.change.emit();

histogram_source.change.emit();
