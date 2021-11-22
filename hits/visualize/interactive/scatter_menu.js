var axes, axis, axis_name, i, j, len, len1, menu, menus, name, ref, ref1, squeeze, suffix;

squeeze = function(possibly_array) {
  var squeezed;
  if (Array.isArray(possibly_array)) {
    squeezed = possibly_array[0];
  } else {
    squeezed = possibly_array;
  }
  return squeezed;
};

menus = {
  'x': x_menu,
  'y': y_menu
};

axes = {
  'x': x_axis,
  'y': y_axis
};

ref = ['x', 'y'];
for (i = 0, len = ref.length; i < len; i++) {
  axis_name = ref[i];
  menu = menus[axis_name];
  axis = axes[axis_name];
  name = squeeze(menu.value);
  scatter_source.data[axis_name] = scatter_source.data[name];
  filtered_source.data[axis_name] = filtered_source.data[name];
  ref1 = ['_all', '_bins_left', '_bins_right'];
  for (j = 0, len1 = ref1.length; j < len1; j++) {
    suffix = ref1[j];
    histogram_source.data[axis_name + suffix] = histogram_source.data[name + suffix];
  }
  axis[0].axis_label = name;
}

// Note: need to re-add calculation of histograms.
scatter_source.change.emit();

filtered_source.change.emit();

histogram_source.change.emit();
