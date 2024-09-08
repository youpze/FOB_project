import csv

with open('HGVS_2020_big_benchmark.tsv') as file:
    dataset = csv.reader(file, delimiter='\t')

    ids = []
    unique_ids = []
    pathogenic = 0
    benign = 0

    for line in dataset:
        if 'Pathogenic' in line:
            pathogenic += 1
        else:
            if 'Benign' in line:
                benign += 1

        for i in line:
            if i.startswith('NC'):
                ids.append(i)

    for id in ids:
        unique_id = id.partition('.')[0]
        if unique_id not in unique_ids:
            unique_ids.append(unique_id)

    # Number of unique ID's
    print('Number of unique IDs is:', len(unique_ids))

    # Number of pathogenic and benign instances
    print('Number of pathogenic instances:', pathogenic)
    print('Number of benign instances:', benign)
