models = cb_obj.document._all_models_by_name._dict
checkbox_groups = (v for k, v of models when k.startsWith('sub_') or k.startsWith('top_'))
for group in checkbox_groups
    group.active = []
    group.trigger('change')

group.callback.func(full_source, group)
