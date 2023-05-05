var all_matches, column, column_names, i, is_case_sensitive, j, len, matches_from_column, matches_query, matches_relevant_query, possibly_lowercase, queries, search_input, target, targets,
  indexOf = [].indexOf;

column_names = {column_names};

search_input = search_input.value;

is_case_sensitive = case_sensitive.active.length > 0;

if (!is_case_sensitive) {
  search_input = search_input.toLowerCase();
  possibly_lowercase = function(t) {
    return t.toString().toLowerCase();
  };
} else {
  possibly_lowercase = function(t) {
    return t;
  };
}

matches_query = function(target) {
  return function(query) {
    return possibly_lowercase(target).indexOf(query) > -1;
  };
};

if (search_input.indexOf('|') > -1) {
  queries = search_input.split('|');
  matches_relevant_query = function(target) {
    return queries.some(matches_query(target));
  };
} else if (search_input.indexOf('&') > -1) {
  queries = search_input.split('&');
  matches_relevant_query = function(target) {
    return queries.every(matches_query(target));
  };
} else {
  matches_relevant_query = function(target) {
    return matches_query(target)(search_input);
  };
}

all_matches = [];

if (search_input.value !== '') {
  for (j = 0, len = column_names.length; j < len; j++) {
    column = column_names[j];
    targets = scatter_source.data[column];
    matches_from_column = (function() {
      var k, len1, results;
      results = [];
      for (i = k = 0, len1 = targets.length; k < len1; i = ++k) {
        target = targets[i];
        if (matches_relevant_query(target) && indexOf.call(all_matches, i) < 0) {
          results.push(i);
        }
      }
      return results;
    })();
    all_matches.push(...matches_from_column);
  }
}

scatter_source.selected.indices = all_matches;
