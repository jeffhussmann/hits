var data, dqs, i, label, len, models, name, ref, search_string;

models = cb_obj.document._all_models_by_name._dict;

models['clear_constraints'].callback.func(models['clear_constraints'], null, require, exports);

search_string = cb_obj.value;

data = models['labels'][0].source.data;

dqs = (function() {
  var i, len, ref, results;
  ref = data['label'];
  results = [];
  for (i = 0, len = ref.length; i < len; i++) {
    name = ref[i];
    results.push(name.indexOf(search_string) === -1);
  }
  return results;
})();

ref = models['labels'];
for (i = 0, len = ref.length; i < len; i++) {
  label = ref[i];
  label.source.data['dq'] = dqs;
}

models['hover_tool'].callback.func(models['hover_tool'], 'force_redraw', require, exports);
