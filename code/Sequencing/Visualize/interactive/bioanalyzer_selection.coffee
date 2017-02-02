models = cb_obj.document._all_models_by_name._dict
sources = (v for k, v of models when k.startsWith('source'))
checkbox_groups = (v for k, v of models when k.startsWith('sub_'))

selected = (s for s in sources when s.selected['0d'].indices.length > 0)
names = (s.name['source_'.length...-'_plotted'.length] for s in selected)

# Clear the selection so that it isn't seen by later callbacks.
s.selected['0d'].indices = [] for s in sources
delete s.selected['0d'].flag for s in sources
delete s.selected['0d'].glyph = null for s in sources

for group in checkbox_groups
    for label, i in group.labels
        if label in names
            if i not in group.active
                group.active.push i
                group.active.sort()
            else
                group.active = (a for a in group.active when a != i)

    group.trigger('change')

group.callback.func(full_source, group)
