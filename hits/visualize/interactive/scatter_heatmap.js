var axis, heatmap_data, hist_data, i, index, j, k, label_data, len, len1, models, name, num_pairs, num_selected, ref, ref1, scatter_data, suffix, x_name, x_names, y_name, y_names;

models = cb_obj.document._all_models_by_name._dict;

scatter_data = models['scatter_source'].data;

label_data = models['filtered_source'].data;

hist_data = models['histogram_source'].data;

heatmap_data = models['heatmap_source'].data;

num_selected = cb_obj.indices.length;

if (num_selected === 0) {
  // Selection was cleared with ESC. Recover the old selection from the axis
  // labels.
  x_name = models['x_axis'].axis_label;
  y_name = models['y_axis'].axis_label;
  num_pairs = heatmap_data['x_name'].length;
  x_names = heatmap_data['x_name'];
  y_names = heatmap_data['y_name'];
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
  name = heatmap_data[axis + '_name'][index];
  scatter_data[axis] = scatter_data[name];
  label_data[axis] = label_data[name];
  ref1 = ['_all', '_bins_left', '_bins_right'];
  for (k = 0, len1 = ref1.length; k < len1; k++) {
    suffix = ref1[k];
    hist_data[axis + suffix] = hist_data[name + suffix];
  }
  models[axis + '_axis'].axis_label = name;
}

// Call to recompute selection histograms.
models['scatter_selection_callback'].func(models['scatter_source'].selected, 'from_heatmap', require, exports);

models['scatter_source'].change.emit();

models['filtered_source'].change.emit();

models['histogram_source'].change.emit();
