var filtered_data, full_data, i, indices, key, models, values;

models = cb_obj.document._all_models_by_name._dict;

full_data = models['scatter_source'].data;

filtered_data = models['labels_source'].data;

indices = cb_obj.selected['1d'].indices;

for (key in full_data) {{
  values = full_data[key];
  filtered_data[key] = (function() {{
    var j, len, results;
    results = [];
    for (j = 0, len = indices.length; j < len; j++) {{
      i = indices[j];
      results.push(values[i]);
    }}
    return results;
  }})();
}}

if ((models['table'] != null)) {{
  models['table'].trigger('change');
}}

models['labels_source'].trigger('change');

requirejs.config({{
    paths: {{
        d3: "http://d3js.org/d3.v4.min"
    }}
}});

requirejs(["d3"], function(d3) {{
    bins = {bins}
    models = cb_obj.document._all_models_by_name._dict

    idxs = cb_obj.selected['1d']['indices'];
    xs = cb_obj.data['x'];
    ys = cb_obj.data['y'];
    xs_selected = [];
    ys_selected = [];

    for (i = 0, len = idxs.length; i < len; i++) {{
        xs_selected.push(xs[idxs[i]]);
        ys_selected.push(ys[idxs[i]]);
    }}

    h = d3.histogram().domain([bins[0], bins[bins.length - 1]]).thresholds(bins);

    xs_binned = h(xs_selected);
    xs_counts = xs_binned.map(function(b) {{ return b.length; }});

    ys_binned = h(ys_selected);
    ys_counts = ys_binned.map(function(b) {{ return b.length; }});

    models['histogram_source'].data['x_selected'] = xs_counts;
    models['histogram_source'].data['y_selected'] = ys_counts;

    models['histogram_source'].trigger('change');
}});
