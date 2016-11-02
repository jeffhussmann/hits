models = cb_obj.document._all_models_by_name._dict

lines = (v for k, v of models when k.startsWith('line_'))

resolution = cb_obj.labels[cb_obj.active][...-(' resolution'.length)]

if resolution == 'nucleotide'
    fig.x_range.start = fig.x_range.start * 3
    fig.x_range.end = fig.x_range.end * 3
else
    fig.x_range.start = fig.x_range.start / 3
    fig.x_range.end = fig.x_range.end / 3

for line in lines
    name = line.name['line_'.length..]
    circle = models['circle_' + name]
    source = models['source_' + name + '_' + resolution]
    line.data_source.data = source.data
    circle.data_source.data = source.data
