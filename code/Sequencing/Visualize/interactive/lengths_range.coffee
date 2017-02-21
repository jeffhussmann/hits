models = cb_obj.document._all_models_by_name._dict

x_range = models['x_range']
y_range = models['y_range']

y_min = -(y_range.end - y_range.start) * 0.005
y_range.start = y_min if y_range.start < y_min

x_range.start = -1 if x_range.start < 0
