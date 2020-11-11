var axis, hist_data, i, j, label_data, len, len1, models, name, ref, ref1, scatter_data, squeeze, suffix;

models = cb_obj.document._all_models_by_name._dict;

scatter_data = models['scatter_source'].data;

label_data = models['filtered_source'].data;

hist_data = models['histogram_source'].data;

squeeze = function(possibly_array) {
  var squeezed;
  if (Array.isArray(possibly_array)) {
    squeezed = possibly_array[0];
  } else {
    squeezed = possibly_array;
  }
  return squeezed;
};

ref = ['x', 'y'];
for (i = 0, len = ref.length; i < len; i++) {
  axis = ref[i];
  name = squeeze(models[axis + '_menu'].value);
  scatter_data[axis] = scatter_data[name];
  label_data[axis] = label_data[name];
  ref1 = ['_all', '_bins_left', '_bins_right'];
  for (j = 0, len1 = ref1.length; j < len1; j++) {
    suffix = ref1[j];
    hist_data[axis + suffix] = hist_data[name + suffix];
  }
  models[axis + '_axis'].axis_label = name;
}

// Call to recompute selection histograms.
models['scatter_selection_callback'].func(models['scatter_source'].selected, cb_data, require, exports);

models['scatter_source'].change.emit();

models['filtered_source'].change.emit();

models['histogram_source'].change.emit();
