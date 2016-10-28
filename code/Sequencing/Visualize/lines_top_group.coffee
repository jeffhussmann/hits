models = cb_obj.document._all_models_by_name._dict
top_name = cb_obj.name[4..]

checkbox_groups = (v for k, v of models when k.startsWith('sub_'))
for group in checkbox_groups
    if group.name[4..] == top_name
        if cb_obj.active.length > 0
            group.active = (i for c, i in group.labels)
        else
            group.active = []
        group.callback.func(invisible_legend, group)
