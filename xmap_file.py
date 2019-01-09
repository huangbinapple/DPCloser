def read_file(file_name):
    """
    Return:
        list[(query_id, ref_id, '+'/'-', alignment_info)]
    """
    result = []
    with open(file_name) as fin:
        for line in map(str.rstrip, filter(lambda x: x[0] != '#', fin)):
            tokens = line.split('\t')
            id_query, id_ref = tokens[1: 3]
            oritation = tokens[7]
            alignment_info = [tuple(map(int, ele.split(','))) for ele in tokens[-1].strip('()').split(')(')]
            result.append((id_query, id_ref, oritation, alignment_info))
    return result

def _test_read_file():
    input_file = 'test.xmap'
    for ele in read_file(input_file):
        print(ele)


if __name__ == '__main__':
    _test_read_file()

