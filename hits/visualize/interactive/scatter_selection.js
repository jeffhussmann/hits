// Make the histograms of all data slightly darker if nothing is selected. 
var binned_to_counts, get_domain_info, i, key, loaded, ref, update_bins, values;

if (cb_obj.indices.length === 0) {
  hist_x_all.glyph.fill_alpha = 0.2;
  hist_y_all.glyph.fill_alpha = 0.2;
} else {
  hist_x_all.glyph.fill_alpha = 0.1;
  hist_y_all.glyph.fill_alpha = 0.1;
}

ref = scatter_source.data;
for (key in ref) {
  values = ref[key];
  filtered_source.data[key] = (function() {
    var j, len, ref1, results;
    ref1 = cb_obj.indices;
    results = [];
    for (j = 0, len = ref1.length; j < len; j++) {
      i = ref1[j];
      results.push(values[i]);
    }
    return results;
  })();
}

filtered_source.change.emit();

table.change.emit();

get_domain_info = function(name) {
  var bins_left, bins_right, bounds, domain_info;
  bins_left = histogram_source.data[name + '_bins_left'];
  bins_right = histogram_source.data[name + '_bins_right'];
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
  var binned, binner, counts, data, domain_info, j, len, name, ref1;
  ref1 = ['x', 'y'];
  for (j = 0, len = ref1.length; j < len; j++) {
    name = ref1[j];
    domain_info = get_domain_info(name);
    // histogram behavior is a little unintuitive.
    // Through trial and error, appears that domain needs to be given min
    // and max data value, and thresholds needs to be given array that
    // include min but not max
    binner = loaded['d3'].histogram().domain(domain_info.bounds).thresholds(domain_info.bins);
    data = filtered_source.data[name];
    binned = binner(data);
    counts = binned_to_counts(binned);
    histogram_source.data[name + '_selected'] = counts;
  }
  return histogram_source.change.emit();
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
