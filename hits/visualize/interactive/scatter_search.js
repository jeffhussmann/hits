var all_matches, column, column_names, i, is_case_sensitive, j, len, matches_any_query, matches_from_column, possibly_lowercase, queries, t, targets,
  indexOf = [].indexOf;

column_names = {column_names};

is_case_sensitive = case_sensitive.active.length > 0;

if (!is_case_sensitive) {
  queries = search_input.value.toLowerCase().split('|');
  possibly_lowercase = function(t) {
    return t.toString().toLowerCase();
  };
} else {
  queries = search_input.value.split('|');
  possibly_lowercase = function(t) {
    return t;
  };
}

matches_any_query = function(t) {
  var matches_query;
  matches_query = function(query) {
    return possibly_lowercase(t).indexOf(query) > -1;
  };
  return queries.some(matches_query);
};

all_matches = [];

if (search_input.value !== '') {
  for (j = 0, len = column_names.length; j < len; j++) {
    column = column_names[j];
    targets = scatter_source.data[column];
    matches_from_column = (function() {
      var k, len1, results;
      results = [];
      for (i = k = 0, len1 = targets.length; k < len1; i = ++k) {
        t = targets[i];
        if (matches_any_query(t) && indexOf.call(all_matches, i) < 0) {
          results.push(i);
        }
      }
      return results;
    })();
    all_matches.push(...matches_from_column);
  }
}

scatter_source.selected.indices = all_matches;
