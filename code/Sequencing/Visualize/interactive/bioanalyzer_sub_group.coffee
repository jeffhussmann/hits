models = cb_obj.document._all_models_by_name._dict

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

if active_names.length == 0
    for line in lines
        name = line.name['line_'.length..]
        line.glyph.line_width = 1
        line.glyph.line_alpha = 0.6
        line.glyph.line_color = "black"
else
    for line in lines
        name = line.name['line_'.length..]

        if name in active_names
            line.glyph.line_color = line.nonselection_glyph.line_color
            line.glyph.line_width = 2
            line.glyph.line_alpha = 0.95
        else
            line.glyph.line_color = "black"
            line.glyph.line_width = 1
            line.glyph.line_alpha = line.nonselection_glyph.line_alpha
            
legend = models['legend']
items = (item for item in legend.all_items when item.label.value in active_names)
items.sort (a, b) -> a.label.value.localeCompare b.label.value
legend.items = items

filtered_data = models['table_source'].data

indices = (i for n, i in full_source.data['sample_name'] when n in active_names)
if not indices?
    indices = []
    
for key of filtered_data
    if key == 'index'
        filtered_data.index = indices
    else
        values = full_source.data[key]
        filtered_data[key] = (values[i] for i in indices)

models['table'].trigger('change')
