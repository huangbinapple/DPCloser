import sys
import getopt
import itertools

import fastg_file
from Bio import SeqIO

def parse_node_long_name(long_name):
    short_name, _ = fastg_file.read_long_name(long_name)
    uid, length, _, is_reverse = fastg_file.read_short_name(short_name)
    uid += ('r' if is_reverse else '')
    return uid, length

class Alignment:

    VALID_THRESHOLD = 0.95

    def __init__(self, alignment_line):
        self.line = alignment_line
        tokens = alignment_line.rstrip().split('\t')
        self.query_id = tokens[0]
        self.subject_id = tokens[1]
        self.query_node_id, query_node_len = parse_node_long_name(
            self.query_id
        )
        identity = float(tokens[2]) / 100
        alignment_length = int(tokens[3])
        # num_error = int(tokens[4])
        self.gap_open = int(tokens[5])
        q_start, q_end, s_start, s_end = map(int, tokens[6:10])
        self.num_delete = alignment_length - (abs(q_start - q_end) + 1)
        self.num_insert = alignment_length - (abs(s_start - s_end) + 1)
        self.e_value = float(tokens[10])
        self.bit_score = float(tokens[11])
        self.start = s_start
        self.end = s_end
        self.left = min(s_start, s_end)
        self.right = max(s_start, s_end)
        self.identity = alignment_length * identity / query_node_len
        self.is_valid = True if\
            self.identity > Alignment.VALID_THRESHOLD else False
        self.is_forward = s_end > s_start

    def adjacent_before(self, alignment, overlap, debug=False):
        """This mean self is adjacent to `alignment` and 
        `self` is on the upper stream of `alignment`."""
        if debug:
            print(self.line.rstrip())
            print(alignment.line.rstrip())
        min_insert = min(self.num_insert, alignment.num_insert)
        min_delete = min(self.num_delete, alignment.num_delete)
        # To be True, two alignment must have the same direction.
        is_same_direction = (self.is_forward == alignment.is_forward)
        shift = overlap - 1
        if not is_same_direction:
            return False
        if self.is_forward:
            if debug:
                print('forward')
                print(self.end)
                print(shift)
                print(alignment.start, min_insert, min_delete)
            return alignment.start - min_insert <= self.end - shift <= \
                alignment.start + min_delete
        else:
            if debug:
                print('backward')
                print(self.end)
                print(shift)
                print(alignment.start, min_delete, min_insert)
            return alignment.start - min_delete <= self.end + shift <= \
                alignment.start + min_insert

def read_file(file_name):
    node_id2alignments = {}
    with open(file_name) as fin:
        for line in filter(lambda x: not x.startswith('#'), fin):
            # Parse a line.
            alignemnt = Alignment(line)
            if not alignemnt.is_valid:
                continue
            if alignemnt.query_node_id in node_id2alignments:
                node_id2alignments[alignemnt.query_node_id].append(
                    alignemnt
                )
            else:
                node_id2alignments[alignemnt.query_node_id] =\
                    [alignemnt]
    return node_id2alignments

def is_adjacent(node_a, node_b, node_id2alignments, overlap,
        debug=False):
    for alignment_a, alignment_b in itertools.product(
            node_id2alignments[node_a], node_id2alignments[node_b]):
        if alignment_a.adjacent_before(alignment_b, overlap, debug):
            return True
    return False

def write_file(output_file, node_id2alignments, nodes, overlap):
    debug = False
    fout = open(output_file, 'w')
    for node_id, node in nodes.items():
        if debug:
            for alignement in node_id2alignments[node_id]:
                fout.write(alignement.line)
        fout.write('NODE: ' + node_id)
        fout.write('\n')
        for child_node, overlap in node.children:
            if debug:
                for alignement in node_id2alignments[child_node.uid]:
                    fout.write(alignement.line)
            fout.write('\t'.join(('', child_node.uid, str(overlap),
                'O' if is_adjacent(node_id, child_node.uid,
                    node_id2alignments, overlap, debug) else 'X')
                ))
            fout.write('\n')
        fout.write('\n')
        debug = False
    fout.close()


def printHelpMessage():
    body = '[-h] <-l overlap len> <fastg file> <blast result> <output>'
    print('python3 {} {}'.format(__file__, body))

def main():
    fastg_file_name = ''
    blast_result_file = ''
    output_file = ''
    overlap_len = None
    options, args = getopt.getopt(sys.argv[1:], 'hl:')
    for option, value in options:
        if option == '-l':
            overlap_len = int(value)
        elif option == '-h':
            printHelpMessage()
            sys.exit()
        else:
            printHelpMessage()
            sys.exit()
    fastg_file_name, blast_result_file, output_file = args

    nodes = fastg_file.build_assembly_graph(fastg_file_name, overlap_len)

    write_file(output_file, read_file(blast_result_file),
        nodes, overlap_len)

if __name__ == '__main__':
    main()
