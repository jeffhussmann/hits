var axis_name, buffer, clean_values, clean_xs, clean_ys, hist_counts, hist_max, i, identical_bins, index_is_clean, j, k, len, log_scale, max, min, number_is_clean, ranges, ref, ref1, xs, ys;

log_scale = {log_scale};

identical_bins = {identical_bins};

ranges = {
  'x': x_range,
  'y': y_range,
  'hist_x': hist_x_range,
  'hist_y': hist_y_range
};

xs = scatter_source.data.x;

ys = scatter_source.data.y;

clean_xs = [];

clean_ys = [];

if (log_scale) {
  number_is_clean = function(x) {
    return x !== 0 && !Number.isNaN(x) && x !== 'NaN';
  };
} else {
  number_is_clean = function(x) {
    return !Number.isNaN(x) && x !== 'NaN';
  };
}

index_is_clean = function(i) {
  return number_is_clean(xs[i]) && number_is_clean(ys[i]);
};

for (i = j = 0, ref = xs.length; (0 <= ref ? j < ref : j > ref); i = 0 <= ref ? ++j : --j) {
  if (index_is_clean(i)) {
    clean_xs.push(xs[i]);
    clean_ys.push(ys[i]);
  }
}

clean_values = {
  'x': clean_xs,
  'y': clean_ys
};

ref1 = ['x', 'y'];
for (k = 0, len = ref1.length; k < len; k++) {
  axis_name = ref1[k];
  min = Math.min(...clean_values[axis_name]);
  max = Math.max(...clean_values[axis_name]);
  if (log_scale) {
    ranges[axis_name].start = min / 2;
    ranges[axis_name].end = max * 2;
  } else {
    buffer = (max - min) * 0.05;
    ranges[axis_name].start = min - buffer;
    ranges[axis_name].end = max + buffer;
  }
  // If identical_bins, want ranges of each histogram to be
  // the same so that heights can be meaningfully compared.
  // If not, maximize dynamic range. 
  if (!identical_bins) {
    hist_counts = histogram_source.data[axis_name + '_all'];
    hist_max = Math.max(...hist_counts);
    ranges['hist_' + axis_name].end = hist_max;
  }
}
