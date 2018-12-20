def read_long_name(long_name):
    """
    Examples of `long_name` format:
        EDGE_1_length_133_cov_25.6667:EDGE_138_length_88469_cov_13.404',EDGE_164_length_174228_cov_14.0904;
        EDGE_1_length_133_cov_25.6667;
    Return:
        tuple: (short_name(str), short_name(tuple))
    """
    # Strip ';'.
    long_name = long_name.rstrip(';')
    # Get node name and children' node names.
    tokens = long_name.split(':')
    node_name = tokens[0]
    if len(tokens) > 1:
        children_names = tuple(tokens[1].split(','))
    else:
        children_names = ()
    # Return result.
    return node_name, children_names

def _test_read_long_name():
    inputs = ("EDGE_1_length_133_cov_25.6667:EDGE_138_length_88469_cov_13.404',EDGE_164_length_174228_cov_14.0904",
              "EDGE_1_length_133_cov_25.6667;")
    for input_ in inputs:
        print(read_long_name(input_))

def read_short_name(short_name):
    """
    Example of 'short_name` format:
        EDGE_1_length_133_cov_25.6667
        EDGE_1_length_133_cov_25.6667'
    Return:
        tuple: (uid(str), length(int), cov(float), is_reversed(bool))
    """
    tokens = short_name.split('_')
    uid = tokens[1]
    length = int(tokens[3])
    if tokens[5][-1] == "'":
        is_reversed = True
        coverage = float(tokens[5][:-1])
    else:
        is_reversed = False
        coverage = float(tokens[5])
    # Return result.
    return (uid, length, coverage, is_reversed)

def _test_read_short_name():
    inputs = ('EDGE_1_length_133_cov_25.6667', "EDGE_1_length_133_cov_25.6667'", "EDGE_1_length_133_cov_25.6667")
    for input_ in inputs:
        print(read_short_name(input_))
 
