models = cb_obj.document._all_models_by_name._dict

# cb_data is used to flag if this callback was triggered manually by
# the search button or subset menu callback's. If not, reset the values of those
# widgets.

if (models['search']?) and cb_data != 'from_search'
    models['search'].value = ''

if (models['subset_menu']?) and cb_data != 'from_subset'
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

get_bins = (name) ->
    bins_left = models['histogram_source'].data[name + '_bins_left']
    bins_right = models['histogram_source'].data[name + '_bins_right']
    # slice() to make copy
    bins = bins_left.slice()
    bins.push(bins_right[bins_right.length - 1])
    return bins

binned_to_counts = (binned) -> (b.length for b in binned)

loaded =
    'd3': if (d3?) then d3 else null

update_bins = () ->
    for name in ['x', 'y']
        bins = get_bins(name)
        bounds = [bins[0], bins[bins.length - 1]]
        binner = loaded['d3'].histogram().domain(bounds).thresholds(bins)
        data = filtered_data[name]
        binned = binner(data)
        counts = binned_to_counts(binned)
        models['histogram_source'].data[name + '_selected'] = counts
    
    models['histogram_source'].trigger('change')

if (d3?)
    update_bins()
    return

if not (window.requirejs?)
    return

`
requirejs.config({{
    paths: {{
        d3: "https://d3js.org/d3-array.v1.min"
    }}
}});

requirejs(['d3'], function(d3) {{
    loaded['d3'] = d3;
    update_bins();
}});
`
