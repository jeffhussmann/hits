# Make the histograms of all data slightly darker if nothing is selected. 
if cb_obj.indices.length == 0
    hist_x_all.glyph.fill_alpha = 0.2
    hist_y_all.glyph.fill_alpha = 0.2
else
    hist_x_all.glyph.fill_alpha = 0.1
    hist_y_all.glyph.fill_alpha = 0.1

for key, values of scatter_source.data
    filtered_source.data[key] = (values[i] for i in cb_obj.indices)

filtered_source.change.emit()

table.change.emit()

get_domain_info = (name) ->
    bins_left = histogram_source.data[name + '_bins_left']
    bins_right = histogram_source.data[name + '_bins_right']
    bounds = [bins_left[0], bins_right[bins_right.length - 1]]
    domain_info =
        bins: bins_left
        bounds: bounds
    return domain_info

binned_to_counts = (binned) -> (b.length for b in binned)

loaded =
    'd3': if d3? then d3 else null

update_bins = () ->
    for name in ['x', 'y']
        domain_info = get_domain_info(name)
        # histogram behavior is a little unintuitive.
        # Through trial and error, appears that domain needs to be given min
        # and max data value, and thresholds needs to be given array that
        # include min but not max
        binner = loaded['d3'].histogram().domain(domain_info.bounds).thresholds(domain_info.bins)
        data = filtered_source.data[name]
        binned = binner(data)
        counts = binned_to_counts(binned)
        histogram_source.data[name + '_selected'] = counts
    
    histogram_source.change.emit()

if d3? and d3.histogram
    update_bins()
    return

if not (window.requirejs?)
    return

`
requirejs.config({
    paths: {
        d3: "https://d3js.org/d3-array.v1.min"
    }
});

requirejs(['d3'], function(d3) {
    loaded['d3'] = d3;
    update_bins();
});
`
