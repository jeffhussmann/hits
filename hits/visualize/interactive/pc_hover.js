var alpha, change, circle_sizes, circle_source, circle_sources, colors, data, dq, dqs, highlight, i, j, k, l, label, len, len1, len2, line_alphas, line_widths, models, ref, text_alphas, v, y;

models = cb_obj.document._all_models_by_name._dict;

data = models['labels'][0].source.data;

if (cb_data === 'force_redraw') {
  change = true;
  highlight = false;
} else {
  y = parseInt(cb_data.renderer.name.slice('line_'.length));
  highlight = cb_data.index.line_indices.length > 0;
  change = data['highlight'][y] !== highlight;
}

if (change) {
  data['highlight'][y] = highlight;
  dqs = data['dq'];
  line_widths = (function() {
    var i, len, results;
    results = [];
    for (i = 0, len = dqs.length; i < len; i++) {
      dq = dqs[i];
      results.push(dq ? 2 : 3);
    }
    return results;
  })();
  circle_sizes = (function() {
    var i, len, results;
    results = [];
    for (i = 0, len = dqs.length; i < len; i++) {
      dq = dqs[i];
      results.push(dq ? 4 : 6);
    }
    return results;
  })();
  colors = (function() {
    var i, len, results;
    results = [];
    for (i = 0, len = dqs.length; i < len; i++) {
      dq = dqs[i];
      results.push('black');
    }
    return results;
  })();
  if (highlight) {
    colors[y] = data['hover_color'][y];
    line_alphas = (function() {
      var i, len, results;
      results = [];
      for (i = 0, len = dqs.length; i < len; i++) {
        dq = dqs[i];
        results.push(dq ? 0.1 : 0.3);
      }
      return results;
    })();
    line_alphas[y] = 0.9;
    text_alphas = (function() {
      var i, len, results;
      results = [];
      for (i = 0, len = dqs.length; i < len; i++) {
        dq = dqs[i];
        results.push(dq ? 0.01 : 0.05);
      }
      return results;
    })();
    text_alphas[y] = 0.95;
    line_widths[y] = 6;
    circle_sizes[y] = 8;
  } else {
    line_alphas = (function() {
      var i, len, results;
      results = [];
      for (i = 0, len = dqs.length; i < len; i++) {
        dq = dqs[i];
        results.push(dq ? 0.1 : 0.7);
      }
      return results;
    })();
    text_alphas = (function() {
      var i, len, results;
      results = [];
      for (i = 0, len = dqs.length; i < len; i++) {
        dq = dqs[i];
        results.push(dq ? 0.05 : 0.9);
      }
      return results;
    })();
  }
  ref = models['labels'];
  for (i = 0, len = ref.length; i < len; i++) {
    label = ref[i];
    label.source.data['text_alpha'] = text_alphas;
    label.source.data['text_color'] = colors;
    label.source.change.emit();
  }
  circle_sources = (function() {
    var results;
    results = [];
    for (k in models) {
      v = models[k];
      if (k.startsWith('source_by_x_')) {
        results.push(v);
      }
    }
    return results;
  })();
  for (j = 0, len1 = circle_sources.length; j < len1; j++) {
    circle_source = circle_sources[j];
    circle_source.data['size'] = circle_sizes;
    circle_source.data['alpha'] = line_alphas;
    circle_source.change.emit();
  }
  for (y = l = 0, len2 = line_alphas.length; l < len2; y = ++l) {
    alpha = line_alphas[y];
    models['line_' + y].glyph.line_alpha = alpha;
    models['line_' + y].glyph.line_width = line_widths[y];
    models['line_' + y].change.emit();
  }
}
