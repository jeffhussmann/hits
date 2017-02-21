models = cb_obj.document._all_models_by_name._dict

alpha = cb_obj.value

flatten = (possibly_arrays) ->
    flat = []
    for possibly_array in possibly_arrays
        if not Array.isArray(possibly_array)
            possibly_array = [possibly_array]
        flat.push possibly_array...
    return flat

line_groups = (v for k, v of models when k.startsWith('line_'))
lines = flatten(line_groups)

checkbox_groups = (v for k, v of models when k.startsWith('sub_'))

active_names = []
for group in checkbox_groups
    labels = (group.labels[i] for i in group.active)
    active_names.push labels...

for line in lines
    for line in lines
        # For every line, put the alpha value where metacodon_sub_group will
        # look for it the next time selection changes.
        line.nonselection_glyph.line_alpha = alpha

        # For non-selected lines, change the current alpha.
        name = line.name['line_'.length..]
        if name not in active_names
            line.glyph.line_alpha = alpha
