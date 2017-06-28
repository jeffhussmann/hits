models = cb_obj.document._all_models_by_name._dict

assignments = models['assignment_menu'].options
assignment = models['assignment_menu'].value
landmark_pair = models['landmark_menu'].value

x_range = models['x_range_left']

lines = (v for k, v of models when k.startsWith('line_'))

maxes = []

if assignments.length == 2 and 'codon' in assignments and 'nucleotide' in assignments
    # Record the original values because metacodon_range.coffee might change them.
    original_start = x_range.start
    original_end = x_range.end

    if assignment == 'nucleotide'
        x_range.start = original_start * 3
        x_range.end = original_end * 3
    else
        x_range.start = original_start / 3
        x_range.end = original_end / 3

for line_group in lines
    if not Array.isArray(line_group)
        line_group = [line_group]
    for line in line_group
        exp_name = line.name['line_'.length..]
        source_name = 'source_' + exp_name + '_' + assignment + '_' + landmark_pair
        source = models[source_name]
        line.data_source.data = source.data
        
        circle_group = models['circle_' + exp_name]

        if (circle_group?)
            if not Array.isArray(circle_group)
                circle_group = [circle_group]
            for circle in circle_group
                circle.data_source.data = source.data

    #y_max = Math.max(source.data['y']...)
    #maxes.push y_max

#overall_max = Math.max(maxes...)
#models['y_range'].end = overall_max * 1.05
