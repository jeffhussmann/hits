models = cb_obj.document._all_models_by_name._dict

# cb_data is used to flag if this callback was triggered manually by
# the search button or subset menu callback's. If not, reset the values of those
# widgets.
if cb_data != 'from_search'
    models['search'].value = ''

if cb_data != 'from_subset'
    models['subset_menu'].value = ''

indices = cb_obj.selected['1d'].indices

# Make the histograms of all data slightly darker if nothing is selected. 
if indices.length == 0
    models['hist_x_all'].glyph.fill_alpha = 0.2
    models['hist_y_all'].glyph.fill_alpha = 0.2
else
    models['hist_x_all'].glyph.fill_alpha = 0.1
    models['hist_y_all'].glyph.fill_alpha = 0.1

full_data = models['scatter_source'].data
filtered_data = models['labels_source'].data

for key, values of full_data
    filtered_data[key] = (values[i] for i in indices)

if (models['table']?)
    models['table'].trigger('change')

models['labels_source'].trigger('change')

bins = {bins}
bounds = [bins[0], bins[bins.length - 1]]

binned_to_counts = (binned) -> (b.length for b in binned)
xs = filtered_data['x']
ys = filtered_data['y']

`
requirejs.config({{
    paths: {{
        d3: "http://d3js.org/d3.v4.min"
    }}
}});

requirejs(["d3"], function(d3) {{
    var histogram, data_to_counts, binned, xs_counts, ys_counts

    histogram = d3.histogram().domain(bounds).thresholds(bins);

    data_to_counts = function(data) {{
        binned = histogram(data);
        return binned_to_counts(binned);
    }};

    xs_counts = data_to_counts(xs);
    ys_counts = data_to_counts(ys);

    models['histogram_source'].data['x_selected'] = xs_counts;
    models['histogram_source'].data['y_selected'] = ys_counts;
    models['histogram_source'].trigger('change');
}});
`
