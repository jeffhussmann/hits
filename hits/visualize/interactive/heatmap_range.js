var font_size, k, l, size;

console.log(labels);

size = cb_obj.end - cb_obj.start;

font_size = Math.min(250, Math.round(3000 / size)).toString() + '%';

for (k in labels) {
  l = labels[k];
  l.text_font_size = font_size;
}

if (cb_obj.start < {lower_bound}) {
  cb_obj.start = {lower_bound};
}

if (cb_obj.end > {upper_bound}) {
  cb_obj.end = {upper_bound};
}
