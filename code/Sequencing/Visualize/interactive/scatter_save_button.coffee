models = cb_obj.document._all_models_by_name._dict
columns = {columns}

filtered_data = models['labels_source'].data

# learned from http://stackoverflow.com/questions/14964035/how-to-export-javascript-array-info-to-csv-on-client-side
csv_content = "data:text/csv;charset=utf-8,"

csv_content += columns.join(',') + '\n'

for i in [0...filtered_data[columns[0]].length]
    line = (filtered_data[column][i] for column in columns).join(',') + '\n'
    csv_content += line

encoded = encodeURI(csv_content)
window.open(encoded)
