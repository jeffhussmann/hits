var binned_to_counts, filtered_data, full_data, get_domain_info, i, indices, key, loaded, models, update_bins, values;

models = cb_obj.document._all_models_by_name._dict;

// cb_data is used to flag if this callback was triggered manually by
// the search button or subset menu callback's. If not, reset the values of those
// widgets.

// To prevent this from erasing a selection that was just made, store indices
// and re-assign them afterwards.
indices = cb_obj.indices;

if (cb_data === 'from_heatmap') {

} else {
  if ((models['search'] != null) && cb_data !== 'from_search') {
    models['search'].value = '';
  }
  if ((models['subset_menu'] != null) && cb_data !== 'from_subset') {
    models['subset_menu'].value = '';
  }
}

cb_obj.indices = indices;

// Make the histograms of all data slightly darker if nothing is selected. 
if (indices.length === 0) {
  models['hist_x_all'].glyph.fill_alpha = 0.2;
  models['hist_y_all'].glyph.fill_alpha = 0.2;
} else {
  models['hist_x_all'].glyph.fill_alpha = 0.1;
  models['hist_y_all'].glyph.fill_alpha = 0.1;
}

full_data = models['scatter_source'].data;

filtered_data = models['filtered_source'].data;

for (key in full_data) {
  values = full_data[key];
  filtered_data[key] = (function() {
    var j, len, results;
    results = [];
    for (j = 0, len = indices.length; j < len; j++) {
      i = indices[j];
      results.push(values[i]);
    }
    return results;
  })();
}

models['filtered_source'].change.emit();

if ((models['table'] != null)) {
  models['table'].change.emit();
}

get_domain_info = function(name) {
  var bins_left, bins_right, bounds, domain_info;
  bins_left = models['histogram_source'].data[name + '_bins_left'];
  bins_right = models['histogram_source'].data[name + '_bins_right'];
  bounds = [bins_left[0], bins_right[bins_right.length - 1]];
  domain_info = {
    bins: bins_left,
    bounds: bounds
  };
  return domain_info;
};

binned_to_counts = function(binned) {
  var b, j, len, results;
  results = [];
  for (j = 0, len = binned.length; j < len; j++) {
    b = binned[j];
    results.push(b.length);
  }
  return results;
};

loaded = {
  'd3': typeof d3 !== "undefined" && d3 !== null ? d3 : null
};

update_bins = function() {
  var binned, binner, counts, data, domain_info, j, len, name, ref;
  ref = ['x', 'y'];
  for (j = 0, len = ref.length; j < len; j++) {
    name = ref[j];
    domain_info = get_domain_info(name);
    // histogram behavior is a little unintuitive.
    // Through trial and error, appears that domain needs to be given min
    // and max data value, and thresholds needs to be given array that
    // include min but not max
    binner = loaded['d3'].histogram().domain(domain_info.bounds).thresholds(domain_info.bins);
    data = filtered_data[name];
    binned = binner(data);
    counts = binned_to_counts(binned);
    models['histogram_source'].data[name + '_selected'] = counts;
  }
  return models['histogram_source'].change.emit();
};

if ((typeof d3 !== "undefined" && d3 !== null) && d3.histogram) {
  update_bins();
  return;
}

if (!(window.requirejs != null)) {
  return;
}


requirejs.config({
    paths: {
        d3: "https://d3js.org/d3-array.v1.min"
    }
});

requirejs(['d3'], function(d3) {
    loaded['d3'] = d3;
    update_bins();
});
;
