import csv

with open('HGVS_2020_big_benchmark.tsv') as file:
    dataset = csv.reader(file, delimiter='\t')

    ids = []
    unique_ids = []

    for line in dataset:
        for i in line:
            if i.startswith('NC'):
                ids.append(i)

    for id in ids:
        unique_id = id.partition('.')[0]
        if unique_id not in unique_ids:
            unique_ids.append(unique_id)

    # Number of unique ID's
    print(len(unique_ids))


