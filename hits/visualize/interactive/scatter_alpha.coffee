value = cb_obj.value
if typeof value == 'string'
    value = parseFloat(value)

scatter.glyph.fill_alpha = value