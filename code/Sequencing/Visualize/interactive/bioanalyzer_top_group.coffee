models = cb_obj.document._all_models_by_name._dict
top_name = cb_obj.name['top_'.length..]

checkbox_groups = (v for k, v of models when k.startsWith('sub_'))
for group in checkbox_groups
    if group.name['top_'.length..] == top_name
        if cb_obj.active.length > 0
            group.active = (i for c, i in group.labels)
        else
            group.active = []

group.callback.func(full_source, group)
