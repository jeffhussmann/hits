models = cb_obj.document._all_models_by_name._dict

lines = (v for k, v of models when k.startsWith('line_'))

maxes = []

for line in lines
    name = line.name['line_'.length..]
    circle = models['circle_' + name]
    source = models['source_' + name + '_' + cb_obj.value]
    line.data_source.data = source.data

    if (circle?)
        circle.data_source.data = source.data

    y_max = Math.max(source.data['y']...)
    maxes.push y_max

overall_max = Math.max(maxes...)
models['y_range'].end = overall_max * 1.05
