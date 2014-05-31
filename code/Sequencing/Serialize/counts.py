from collections import Counter

def write_file(counter, file_name):
    with open(file_name, 'w') as fh:
        for string, count in counter.most_common():
            fh.write('{0}\t{1}\n'.format(string, count))

def read_file(file_name):
    counter = Counter()
    for line in open(file_name):
        string, count = line.strip().split()
        counter[string] = int(count)
    return counter

def combine_data(first_counter, second_counter):
    first_counter.update(second_counter)
    return first_counter
