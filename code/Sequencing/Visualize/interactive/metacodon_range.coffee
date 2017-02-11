models = cb_obj.document._all_models_by_name._dict

chooser = models['highest_level_chooser']
resolution = chooser.labels[chooser.active][...-(' resolution'.length)]
if resolution == 'nucleotide'
    x_max = 300
else
    x_max = 100

x_range = models['x_range']
y_range = models['y_range']

y_range.start = 0 if y_range.start < 0
y_range.end = 50 if y_range.end > 50

x_range.start = -x_max if x_range.start < -x_max
x_range.end = x_max if x_range.end > x_max
