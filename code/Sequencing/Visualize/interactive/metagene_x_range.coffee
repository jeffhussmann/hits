models = cb_obj.document._all_models_by_name._dict

resolution_model = models['resolution']
resolution = resolution_model.labels[resolution_model.active][...-(' resolution'.length)]

if resolution == 'nucleotide'
    max_before = 270
    max_after = 750
else
    max_before = 90
    max_after = 250

changed = cb_obj
start_changed = changed.attributes.start != changed._previousAttributes.start
end_changed = changed.attributes.end != changed._previousAttributes.end

changed_landmark = changed.name['x_range_'.length..]
other_landmark = if changed_landmark == 'stop_codon' then 'start_codon' else 'stop_codon'
other = models['x_range_' + other_landmark]

if changed_landmark == 'start_codon'
    changed.start = -max_before if changed.start < -max_before
    changed.end = max_after if changed.end > max_after
else
    changed.end = max_before if changed.end > max_before
    changed.start = -max_after if changed.start < -max_after

if end_changed
    other.start = -changed.end
if start_changed
    other.end = -changed.start
