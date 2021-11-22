size = cb_obj.end - cb_obj.start
font_size = Math.min(250, Math.round(3000 / size)).toString() + '%'
l.text_font_size = font_size for k, l of labels
cb_obj.start = {lower_bound} if cb_obj.start < {lower_bound}
cb_obj.end = {upper_bound} if cb_obj.end > {upper_bound}