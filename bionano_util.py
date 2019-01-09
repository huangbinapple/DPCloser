def read_key_file(file_name):
    result = {}
    with open(file_name) as f:
        for _ in range(4):
            f.readline()
        for line in f:
            tokens = line.split('\t')
            result[tokens[0]] = tokens[1]
    return result