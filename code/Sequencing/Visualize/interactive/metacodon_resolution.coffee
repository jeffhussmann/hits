models = cb_obj.document._all_models_by_name._dict

lines = (v for k, v of models when k.startsWith('line_'))

resolution = cb_obj.value

x_range = models['x_range_left']

console.log 'in resolution'

# Record the original values because metacodon_range.coffee might change them.
original_start = x_range.start
original_end = x_range.end

if resolution == 'nucleotide'
    x_range.start = original_start * 3
    x_range.end = original_end * 3
else
    x_range.start = original_start / 3
    x_range.end = original_end / 3

for line_group in lines
    if not Array.isArray(line_group)
        line_group = [line_group]
    for line in line_group
        name = line.name['line_'.length..]
        source = models['source_' + name + '_' + resolution]
        line.data_source.data = source.data

if (models['x_axis']?)
    models['x_axis'].axis_label = 'Offset (' + resolution + 's)'
