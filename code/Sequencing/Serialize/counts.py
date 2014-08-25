from collections import Counter

def write_file(counter, file_name):
    with open(file_name, 'w') as fh:
        for key, count in counter.most_common():
            if isinstance(key, tuple):
                key = ' '.join(key)
            fh.write('{0}\t{1}\n'.format(key, count))

def read_file(file_name):
    counter = Counter()
    for line in open(file_name):
        key, count = line.strip().split('\t')
        counter[key] = int(count)
    return counter

def combine_data(first_counter, second_counter):
    first_counter.update(second_counter)
    return first_counter
