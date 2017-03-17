var filtered_data, full_data, i, indices, key, models, values, matches_subset;

models = cb_obj.document._all_models_by_name._dict;

full_data = models['scatter_source'].data;

filtered_data = models['labels_source'].data;

indices = cb_obj.selected['1d'].indices;

subset_indices = {subset_indices}

matches_subset = false;

for (key in subset_indices) {{
    if (_.isEqual(indices, subset_indices[key])) {{
        models['subset_menu'].value = key;
        matches_subset = true;
    }}
}}

if (matches_subset == false) {{
    models['subset_menu'].value = '';
}}

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

if (indices.length == 0) {{
    models['hist_x_all'].glyph.fill_alpha = 0.2;
    models['hist_y_all'].glyph.fill_alpha = 0.2;
}} else {{
    models['hist_x_all'].glyph.fill_alpha = 0.1;
    models['hist_y_all'].glyph.fill_alpha = 0.1;
}}

requirejs.config({{
    paths: {{
        d3: "http://d3js.org/d3.v4.min"
    }}
}});

requirejs(["d3"], function(d3) {{
    bins = {bins}

    xs = cb_obj.data['x'];
    ys = cb_obj.data['y'];

    xs_selected = [];
    ys_selected = [];

    for (i = 0, len = indices.length; i < len; i++) {{
        xs_selected.push(xs[indices[i]]);
        ys_selected.push(ys[indices[i]]);
    }}

    bounds = [bins[0], bins[bins.length - 1]];
    h = d3.histogram().domain(bounds).thresholds(bins);

    binned_to_counts = function(binned) {{
        return binned.map(function(b) {{ return b.length; }});
    }};

    xs_binned = h(xs_selected);
    xs_counts = binned_to_counts(xs_binned);

    ys_binned = h(ys_selected);
    ys_counts = binned_to_counts(ys_binned);

    models['histogram_source'].data['x_selected'] = xs_counts;
    models['histogram_source'].data['y_selected'] = ys_counts;

    models['histogram_source'].trigger('change');
}});
