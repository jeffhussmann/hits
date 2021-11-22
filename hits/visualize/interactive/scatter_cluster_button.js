var axis, cluster, data, formatters, i, index, j, k, l, len, len1, num_pairs, order_key, ref, ref1, x_name, x_names, y_name, y_names;

formatters = {
  'clustered': {clustered_formatter},
  'original': {original_formatter}
};

cluster = !dendrogram_lines.visible;

dendrogram_lines.visible = cluster;

if (cluster) {
  order_key = 'clustered';
} else {
  order_key = 'original';
}

data = heatmap_source.data;

ref = ['r', 'x', 'x_name', 'y', 'y_name', 'color'];
for (j = 0, len = ref.length; j < len; j++) {
  k = ref[j];
  data[k] = data[k + '_' + order_key];
}

ref1 = [heatmap_x_axis, heatmap_y_axis];
for (l = 0, len1 = ref1.length; l < len1; l++) {
  axis = ref1[l];
  // Note: axis is a splattable list
  axis[0].formatter.code = formatters[order_key];
}

// Determine the new selection from the scatter axis labels.
// Note: each axis is a splattable list
x_name = x_axis[0].axis_label;

y_name = y_axis[0].axis_label;

num_pairs = data['x_name'].length;

x_names = data['x_name'];

y_names = data['y_name'];

index = ((function() {
  var m, ref2, results;
  results = [];
  for (i = m = 0, ref2 = num_pairs; (0 <= ref2 ? m <= ref2 : m >= ref2); i = 0 <= ref2 ? ++m : --m) {
    if (x_names[i] === x_name && y_names[i] === y_name) {
      results.push(i);
    }
  }
  return results;
})())[0];

heatmap_source.selected.indices = [index];

heatmap_source.change.emit();
