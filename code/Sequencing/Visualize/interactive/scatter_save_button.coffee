models = cb_obj.document._all_models_by_name._dict
columns = {columns}

filtered_data = models['labels_source'].data

# learned from http://stackoverflow.com/questions/14964035/how-to-export-javascript-array-info-to-csv-on-client-side
csv_content = "data:text/csv;charset=utf-8,"

lines = [columns.join('\t')]

for i in [0...filtered_data[columns[0]].length]
    line = (filtered_data[column][i] for column in columns).join('\t')
    lines.push line

csv_content += lines.join('\n')

encoded = encodeURI(csv_content)
link = document.createElement('a')
link.setAttribute('href', encoded)
link.setAttribute('download', 'table_contents.txt')
link.click()
