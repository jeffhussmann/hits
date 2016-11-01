colors_dict = {colors_dict}

models = cb_obj.document._all_models_by_name._dict

lines = (v for k, v of models when k.startsWith('line_'))
circles = (v for k, v of models when k.startsWith('circle_'))
checkbox_groups = (v for k, v of models when k.startsWith('sub_'))

active_names = []
for group in checkbox_groups
    labels = (group.labels[i] for i in group.active)
    active_names.push labels...

if active_names.length == 0
    for line in lines
        name = line.name['line_'.length..]
        color = colors_dict[name]
        line.glyph.line_width = 1
        line.glyph.line_alpha = 0.6
        line.glyph.line_color = "black"
        
    circle.glyph.visible = false for circle in circles
        
else
    for line in lines
        name = line.name['line_'.length..]
        color = colors_dict[name]

        if name in active_names
            line.glyph.line_color = color
            line.hover_glyph.line_color = color
            line.glyph.line_width = 2
            line.glyph.line_alpha = 0.95
        else
            line.glyph.line_color = "black"
            line.glyph.line_width = 1
            line.glyph.line_alpha = 0.2
            
    for circle in circles
        name = circle.name['circle_'.length..]
        color = colors_dict[name]

        if name in active_names
            circle.glyph.visible = true
            circle.glyph.line_color = color
            circle.glyph.fill_color = color
        else
            circle.glyph.visible = false
        
legend = (v for k, v of models when k == 'legend')[0]
legend.items = (item for item in invisible_legend.items when item.label.value in active_names)
