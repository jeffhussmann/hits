var all_matches, column, column_names, i, is_case_sensitive, j, len, matches, possibly_lowercase, query, t, targets,
  indexOf = [].indexOf;

column_names = {column_names};

is_case_sensitive = case_sensitive.active.length > 0;

if (!is_case_sensitive) {
  query = search_input.value.toLowerCase();
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
    targets = scatter_source.data[column];
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

scatter_source.selected.indices = all_matches;
