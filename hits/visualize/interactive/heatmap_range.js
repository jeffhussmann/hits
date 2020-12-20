var font_size, i, l, len, models, ref, size;

models = cb_obj.document._all_models_by_name._dict;

size = cb_obj.end - cb_obj.start;

font_size = Math.min(250, Math.round(3000 / size)).toString() + '%';

ref = models['labels'];
for (i = 0, len = ref.length; i < len; i++) {
  l = ref[i];
  l.text_font_size = font_size;
}

if (cb_obj.start < {lower_bound}) {
  cb_obj.start = {lower_bound};
}

if (cb_obj.end > {upper_bound}) {
  cb_obj.end = {upper_bound};
}
