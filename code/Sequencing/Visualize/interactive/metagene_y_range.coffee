models = cb_obj.document._all_models_by_name._dict

y_range = models['y_range']
y_range.start = 0 if y_range.start < 0
y_range.end = 50 if y_range.end > 50
