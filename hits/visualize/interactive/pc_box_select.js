var dqs, highlight, i, j, k, l, labels, len, len1, len2, line, lines, models, ref, v, y, ys;

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

ys = [];

for (i = 0, len = lines.length; i < len; i++) {
  line = lines[i];
  y = line.name.slice('line_'.length);
  highlight = !line.data_source.data['dq'].some(function(x) {
    return x;
  });
  if (highlight) {
    ys.push(y);
  }
}

dqs = (function() {
  var j, len1, results;
  results = [];
  for (j = 0, len1 = lines.length; j < len1; j++) {
    line = lines[j];
    results.push(true);
  }
  return results;
})();

for (j = 0, len1 = ys.length; j < len1; j++) {
  y = ys[j];
  dqs[y] = false;
}

ref = models['labels'];
for (l = 0, len2 = ref.length; l < len2; l++) {
  labels = ref[l];
  labels.source.data['dq'] = dqs;
}

models['constraints'].change.emit();

models['hover_tool'].callback.func(models['hover_tool'], 'force_redraw', require, exports);
