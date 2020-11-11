var bottom, by_x_dq, by_y_dqs, dq, i, index, indices, j, len, models, name, source_by_y, top, val, vals, vbar_data, x, x_int, y;

models = cb_obj.document._all_models_by_name._dict;

by_y_dqs = (function() {
  var results;
  results = [];
  for (name in models) {
    source_by_y = models[name];
    if (name.startsWith('source_by_y')) {
      results.push(source_by_y.data['dq']);
    }
  }
  return results;
})();

x = cb_obj.name.slice('source_by_x_'.length);

indices = cb_obj.selected.indices;

if (indices.length > 0) {
  by_x_dq = (function() {
    var j, len, ref, results;
    ref = cb_obj.data['dq'];
    results = [];
    for (j = 0, len = ref.length; j < len; j++) {
      val = ref[j];
      results.push(true);
    }
    return results;
  })();
  for (j = 0, len = indices.length; j < len; j++) {
    i = indices[j];
    by_x_dq[i] = false;
  }
  cb_obj.data['dq'] = by_x_dq;
  for (y in by_x_dq) {
    dq = by_x_dq[y];
    by_y_dqs[y][x] = dq;
  }
  // Don't show a bar if all points in the column are in the selection.
  if (indices.length === cb_obj.data['dq'].length) {
    bottom = 0;
    top = 0;
  } else {
    vals = (function() {
      var k, len1, results;
      results = [];
      for (k = 0, len1 = indices.length; k < len1; k++) {
        i = indices[k];
        results.push(cb_obj.data['y'][i]);
      }
      return results;
    })();
    bottom = Math.min(...vals) - 10;
    top = Math.max(...vals) + 10;
  }
  x_int = parseInt(x);
  vbar_data = models['constraints'].data;
  index = vbar_data['x'].indexOf(x_int);
  if (index === -1) {
    index = vbar_data['x'].length;
  }
  vbar_data['x'][index] = x_int;
  vbar_data['top'][index] = top;
  vbar_data['bottom'][index] = bottom;
  vbar_data['width'][index] = 0.2;
}

cb_obj.selected.indices = [];

models['search'].value = '';
