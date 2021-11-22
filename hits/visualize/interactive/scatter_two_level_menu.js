var axis_name, axis_subtitle, axis_title, i, j, len, len1, menu_0, menu_1, name, ref, ref1, subtitle, suffix, title;

ref = ['x', 'y'];
for (i = 0, len = ref.length; i < len; i++) {
  axis_name = ref[i];
  if (axis_name === 'x') {
    menu_0 = x_menu_0;
    menu_1 = x_menu_1;
    axis_title = x_axis_title;
    axis_subtitle = x_axis_subtitle;
  } else {
    menu_0 = y_menu_0;
    menu_1 = y_menu_1;
    axis_title = y_axis_title;
    axis_subtitle = y_axis_subtitle;
  }
  title = menu_0.value[0];
  subtitle = menu_1.value[0];
  name = title + ' ' + subtitle;
  scatter_source.data[axis_name] = scatter_source.data[name];
  filtered_source.data[axis_name] = filtered_source.data[name];
  ref1 = ['_all', '_bins_left', '_bins_right'];
  for (j = 0, len1 = ref1.length; j < len1; j++) {
    suffix = ref1[j];
    histogram_source.data[axis_name + suffix] = histogram_source.data[name + suffix];
  }
  axis_title.text = title;
  axis_subtitle.text = subtitle;
}

// Note: need to re-add calculation of histograms.
scatter_source.change.emit();

filtered_source.change.emit();

histogram_source.change.emit();
