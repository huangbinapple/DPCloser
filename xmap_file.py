def read_file(file_name):
    """
    Return:
        list[(query_id, ref_id, '+'/'-', alignment_info)]
    Notes:
        Indexs in `alignment info` is 1-indexed.
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

def read_file_2(file_name):
    """
    Notes:
        All indexs are 0-indexed.
    """
    alignments = []
    for (id_query, id_ref, oritation, alignment_info) in \
        read_file(file_name):
        alignment = {}
        alignment['id_query'] = id_query
        alignment['id_ref'] = id_ref
        alignment['oritation'] = True if oritation == '+' else False
        alignment['subject_start_index'] = alignment_info[0][0] - 1
        alignment['query_start_index'] = alignment_info[0][1] - 1
        alignment['subject_end_index'] = alignment_info[-1][0] - 1
        alignment['query_end_index'] = alignment_info[-1][1] - 1
        if not alignment['oritation']:
            alignment['query_start_index'], alignment['query_end_index'] =\
                alignment['query_end_index'], alignment['query_start_index'] 
        alignments.append(alignment)
    return alignments

def _test_read_file():
    input_file = 'test.xmap'
    for ele in read_file(input_file):
        print(ele)


if __name__ == '__main__':
    _test_read_file()

