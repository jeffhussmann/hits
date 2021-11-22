rect_data = rect_source.data
length = rect_data['x'].length
xs = rect_data['x']
ys = rect_data['y']
widths = rect_data['width']
heights = rect_data['height']

cross_data = cross_source.data
cross_xs = cross_data['x']
cross_ys = cross_data['y']

# learned from http://stackoverflow.com/questions/14964035/how-to-export-javascript-array-info-to-csv-on-client-side
csv_content = "data:text/csv;charset=utf-8,"

lines = []

for i in [0...length]
    x0 = xs[i] - widths[i] / 2
    x1 = xs[i] + widths[i] / 2
    y0 = ys[i] - heights[i] / 2
    y1 = ys[i] + heights[i] / 2
    line = [x0, x1, y0, y1, cross_xs[i], cross_ys[i]].join(',')
    lines.push line

csv_content += lines.join('\n')

encoded = encodeURI(csv_content)
link = document.createElement('a')
link.setAttribute('href', encoded)
link.setAttribute('download', 'rectangles_{{title}}.txt')
link.click()
