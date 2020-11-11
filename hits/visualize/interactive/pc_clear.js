var c, circle_source, i, k, l, label, len, line_source, lines, models, name, ref, v;

models = cb_obj.document._all_models_by_name._dict;

lines = (function() {
  var results;
  results = [];
  for (k in models) {
    v = models[k];
    if (k.startsWith('line')) {
      results.push(v);
    }
  }
  return results;
})();

for (name in models) {
  line_source = models[name];
  if (name.startsWith('source_by_y')) {
    line_source.data['dq'] = (function() {
      var i, len, ref, results;
      ref = line_source.data['dq'];
      results = [];
      for (i = 0, len = ref.length; i < len; i++) {
        c = ref[i];
        results.push(false);
      }
      return results;
    })();
  }
}

for (name in models) {
  circle_source = models[name];
  if (k.startsWith('source_by_x_')) {
    circle_source.data['dq'] = (function() {
      var i, len, results;
      results = [];
      for (i = 0, len = lines.length; i < len; i++) {
        l = lines[i];
        results.push(false);
      }
      return results;
    })();
  }
}

ref = models['labels'];
for (i = 0, len = ref.length; i < len; i++) {
  label = ref[i];
  label.source.data['dq'] = (function() {
    var j, len1, results;
    results = [];
    for (j = 0, len1 = lines.length; j < len1; j++) {
      l = lines[j];
      results.push(false);
    }
    return results;
  })();
}

models['constraints'].data = {
  'x': [],
  'width': [],
  'top': [],
  'bottom': []
};

models['constraints'].change.emit();

models['hover_tool'].callback.func(models['hover_tool'], 'force_redraw', require, exports);
