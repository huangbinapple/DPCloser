import sys
import getopt
import itertools

import fastg_file
import dot_file
from Bio import SeqIO

def parse_node_long_name(long_name):
    short_name, _ = fastg_file.read_long_name(long_name)
    uid, length, _, is_reverse = fastg_file.read_short_name(short_name)
    uid += ('r' if is_reverse else '')
    return uid, length

class Alignment:

    VALID_THRESHOLD = 0.95
    ERROR_MARGIN = 1
    COLOR_PROFILE = ['green', 'yellow', 'yellow', 'orange', 'orange',
        'red']

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
        self.num_mismatch = int(tokens[4])
        self.gap_open = int(tokens[5])
        q_start, q_end, s_start, s_end = map(int, tokens[6:10])
        self.start_cut = q_start - 1
        self.end_cut = query_node_len - q_end
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
        self.children = []
        self.color = self.COLOR_PROFILE[min(len(self.COLOR_PROFILE) - 1,
            self.num_mismatch + self.gap_open)]

    def __str__(self):
        return ','.join((self.query_node_id,
            '-'.join(map(str, (self.start, self.end))),
            str(self.end - self.start + 1)
            ))

    def add_child(self, alignment):
        self.children.append(alignment)

    def adjacent_before(self, alignment, overlap):
        """This mean self is adjacent to `alignment` and 
        `self` is on the upper stream of `alignment`."""
        min_insert = min(self.num_insert, alignment.num_insert)
        min_delete = min(self.num_delete, alignment.num_delete)
        # To be True, two alignment must be in the same ref and have
        # the same direction.
        is_same_dir = (self.is_forward == alignment.is_forward and
            self.subject_id == alignment.subject_id)
        shift = overlap - 1
        if not is_same_dir:
            return False
        if self.is_forward:
            real_self_end = self.end + self.end_cut
            result_other_start = alignment.start - alignment.start_cut
            valid_l = result_other_start  - min_insert - \
                Alignment.ERROR_MARGIN
            valid_h = result_other_start + min_delete + \
                Alignment.ERROR_MARGIN
            return valid_l <= real_self_end - shift <= valid_h
        else:
            real_self_end = self.end - self.end_cut
            real_other_start = alignment.start + alignment.start_cut
            valid_l = real_other_start - min_delete - \
                Alignment.ERROR_MARGIN
            valid_h = real_other_start + min_insert + \
                Alignment.ERROR_MARGIN
            return valid_l <= real_self_end + shift <= valid_h

    @classmethod
    def index(cls, alignments, key):
        if key == 'node id':
            key_attr = 'query_node_id'
        elif key == 'start position':
            key_attr = 'start'
        else:
            raise ValueError('Key not supported.')
        index = {}
        for alignment in alignments:
            if getattr(alignment, key_attr) in index:
                index[getattr(alignment, key_attr)].append(alignment)
            else:
                index[getattr(alignment, key_attr)] = [alignment]
        return index

    @classmethod
    def add_connection(cls, alignments, nodes):
        """Only for those forward alignments"""
        alignments = list(filter(lambda x: x.is_forward, alignments))
        # Sort method.
        alignments.sort(key=lambda x: x.start)
        for i in range(len(alignments)):
            for child_node, overlap in\
                    nodes[alignments[i].query_node_id].children:
                safe_margin = 50
                stop_position = alignments[i].end - (overlap - 1) + \
                    safe_margin
                j = i + 1
                while j < len(alignments) and \
                        alignments[j].start < stop_position:
                    if alignments[j].query_node_id == \
                            child_node.uid and \
                            alignments[i].adjacent_before(alignments[j],
                                overlap):
                        alignments[i].add_child(alignments[j])
                    j += 1

    @classmethod
    def write_alignments_to_dot_file(cls, alignments, file_name,
            actions=None):
        with dot_file.DotFile(file_name) as fout:
            for alignment in alignments:
                fout.add_node(str(alignment),
                    {"color": alignment.color})
                if alignment.children:
                    for child in alignment.children:
                        attribute = {}
                        if actions and actions[alignment] == child:
                            attribute['color'] = 'green'
                        fout.add_edge(*map(str, (alignment, child)),
                            attribute)

    @classmethod
    def get_path(cls, alignments):
        values = {a: (1 if a.color == 'green' else 0) \
            for a in alignments}
        actions = {a: (a.children[0] if a.children else None) \
            for a in alignments}
        for _ in range(len(alignments)):
            for alignment in alignments:
                if actions[alignment]:
                    # Update values of alignments.
                    values[alignment] = values[actions[alignment]] + \
                        1 if alignment.color == 'green' else 0

                    children_values = [values[child] for \
                        child in alignment.children]
                    # Update action of alignments. 
                    actions[alignment] = alignment.children[
                        children_values.index(max(children_values))
                    ]
        return values, actions

def read_file(file_name):
    alignments = []
    with open(file_name) as fin:
        for line in filter(lambda x: not x.startswith('#'), fin):
            # Parse a line.
            alignment = Alignment(line)
            alignments.append(alignment)
    return alignments

def is_adjacent(node_a, node_b, node_id2alignments, overlap):
    for alignment_a, alignment_b in itertools.product(
            node_id2alignments[node_a], node_id2alignments[node_b]):
        if alignment_a.adjacent_before(alignment_b, overlap):
            return True
    return False

def write_file(output_file, alignments):
    fout = open(output_file, 'w')
    for alignment in alignments:
        fout.write(alignment.line)
        for child in alignment.children:
            fout.write('\t')
            fout.write(child.line)
        fout.write('\n')
    fout.close()

def print_help_message():
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
            print_help_message()
            sys.exit()
        else:
            print_help_message()
            sys.exit()
    fastg_file_name, blast_result_file, output_file = args

    nodes = fastg_file.build_assembly_graph(fastg_file_name,
        overlap_len)
    alignments = list(filter(lambda x: x.is_valid and x.is_forward,
        read_file(blast_result_file)))
    Alignment.add_connection(alignments, nodes)
    alignments.sort(key=lambda x: x.start)
    # write_file(output_file, alignments)
    actions = Alignment.get_path(alignments)[1]
    Alignment.write_alignments_to_dot_file(alignments, output_file,
        actions)

if __name__ == '__main__':
    main()
