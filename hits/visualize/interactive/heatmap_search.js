var c, col_alphas, col_labels, col_matches, data, i, j, k, len, len1, query, r, row_alphas, row_labels, row_matches;

data = quad_source.data;

row_labels = {row_labels};

col_labels = {col_labels};

query = cb_obj.value;

if (query !== '') {
  row_matches = (function() {
    var j, len, results;
    results = [];
    for (i = j = 0, len = row_labels.length; j < len; i = ++j) {
      r = row_labels[i];
      if (r.indexOf(query) > -1) {
        results.push(i);
      }
    }
    return results;
  })();
  row_alphas = (function() {
    var j, len, results;
    results = [];
    for (j = 0, len = row_labels.length; j < len; j++) {
      r = row_labels[j];
      results.push(0.3);
    }
    return results;
  })();
  for (j = 0, len = row_matches.length; j < len; j++) {
    i = row_matches[j];
    row_alphas[i] = 1;
  }
  col_matches = (function() {
    var k, len1, results;
    results = [];
    for (i = k = 0, len1 = col_labels.length; k < len1; i = ++k) {
      c = col_labels[i];
      if (c.indexOf(query) > -1) {
        results.push(i);
      }
    }
    return results;
  })();
  col_alphas = (function() {
    var k, len1, results;
    results = [];
    for (k = 0, len1 = col_labels.length; k < len1; k++) {
      c = col_labels[k];
      results.push(0.3);
    }
    return results;
  })();
  for (k = 0, len1 = col_matches.length; k < len1; k++) {
    i = col_matches[k];
    col_alphas[i] = 1;
  }
} else {
  row_matches = [];
  row_alphas = (function() {
    var l, len2, results;
    results = [];
    for (l = 0, len2 = row_labels.length; l < len2; l++) {
      r = row_labels[l];
      results.push(1);
    }
    return results;
  })();
  col_matches = [];
  col_alphas = (function() {
    var l, len2, results;
    results = [];
    for (l = 0, len2 = col_labels.length; l < len2; l++) {
      c = col_labels[l];
      results.push(1);
    }
    return results;
  })();
}

data['left'] = [
  ...((function() {
    var l,
  len2,
  results;
    results = [];
    for (l = 0, len2 = row_matches.length; l < len2; l++) {
      i = row_matches[l];
      results.push({lower_bound});
    }
    return results;
  })()),
  ...((function() {
    var l,
  len2,
  results;
    results = [];
    for (l = 0, len2 = col_matches.length; l < len2; l++) {
      i = col_matches[l];
      results.push({lower_bound} + i);
    }
    return results;
  })())
];

data['right'] = [
  ...((function() {
    var l,
  len2,
  results;
    results = [];
    for (l = 0, len2 = row_matches.length; l < len2; l++) {
      i = row_matches[l];
      results.push({x_upper_bound});
    }
    return results;
  })()),
  ...((function() {
    var l,
  len2,
  results;
    results = [];
    for (l = 0, len2 = col_matches.length; l < len2; l++) {
      i = col_matches[l];
      results.push({lower_bound} + i + 1);
    }
    return results;
  })())
];

data['bottom'] = [
  ...((function() {
    var l,
  len2,
  results;
    results = [];
    for (l = 0, len2 = row_matches.length; l < len2; l++) {
      i = row_matches[l];
      results.push({y_upper_bound} - i - 1);
    }
    return results;
  })()),
  ...((function() {
    var l,
  len2,
  results;
    results = [];
    for (l = 0, len2 = col_matches.length; l < len2; l++) {
      i = col_matches[l];
      results.push({lower_bound});
    }
    return results;
  })())
];

data['top'] = [
  ...((function() {
    var l,
  len2,
  results;
    results = [];
    for (l = 0, len2 = row_matches.length; l < len2; l++) {
      i = row_matches[l];
      results.push({y_upper_bound} - i);
    }
    return results;
  })()),
  ...((function() {
    var l,
  len2,
  results;
    results = [];
    for (l = 0, len2 = col_matches.length; l < len2; l++) {
      i = col_matches[l];
      results.push({y_upper_bound});
    }
    return results;
  })())
];

row_label_source.data['alpha'] = row_alphas;

col_label_source.data['alpha'] = col_alphas;

quad_source.change.emit();

row_label_source.change.emit();

col_label_source.change.emit();
