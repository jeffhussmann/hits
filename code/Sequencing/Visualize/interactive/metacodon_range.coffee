models = cb_obj.document._all_models_by_name._dict
x_range = models['x_range']
y_range = models['y_range']

y_range.start = 0 if y_range.start < 0
y_range.end = 50 if y_range.end > 50

x_range.start = -100 if x_range.start < -100
x_range.end = 100 if x_range.end > 100
