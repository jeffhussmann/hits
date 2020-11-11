var all_matches, case_sensitive, column, column_names, i, j, len, matches, models, possibly_lowercase, query, t, targets,
  indexOf = [].indexOf;

models = cb_obj.document._all_models_by_name._dict;

column_names = {column_names};

case_sensitive = models['case_sensitive'].active.length > 0;

if (!case_sensitive) {
  query = cb_obj.value.toLowerCase();
  possibly_lowercase = function(t) {
    return t.toString().toLowerCase();
  };
} else {
  query = cb_obj.value;
  possibly_lowercase = function(t) {
    return t;
  };
}

all_matches = [];

if (query !== '') {
  for (j = 0, len = column_names.length; j < len; j++) {
    column = column_names[j];
    targets = models['scatter_source'].data[column];
    matches = (function() {
      var k, len1, results;
      results = [];
      for (i = k = 0, len1 = targets.length; k < len1; i = ++k) {
        t = targets[i];
        if (possibly_lowercase(t).indexOf(query) > -1 && indexOf.call(all_matches, i) < 0) {
          results.push(i);
        }
      }
      return results;
    })();
    all_matches.push(...matches);
  }
}

models['scatter_source'].selected.indices = all_matches;

//models['scatter_selection_callback'].func(models['scatter_source'].selected, 'from_search', require, exports)
