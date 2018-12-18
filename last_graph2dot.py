import sys
import getopt
import subprocess
import math

import dotFile

BASE = 60
K = 0.5
def calc_node_size(seq_len):
    return math.log(seq_len, BASE) + K

def write_dot_file(input_file, output_file):
    # Extract edges.
    edges = []
    seq_lens = {}
    with open(input_file) as fin:
        for line in fin:
            if line.startswith('ARC'):
                tokens = line.split('\t')
                edge = tuple(map(int, tokens[1:3]))
                edge_reverse = (-edge[1], -edge[0])
                edges.extend((edge, edge_reverse))
            elif line.startswith('NODE'):
                tokens = line.split('\t', 3)
                node_id = tokens[1]
                seq_lens[node_id] = int(tokens[2])
            else:
                continue

    # Write dot.
    with dotFile.DotFile(output_file) as fout:
        # Write nodes.
        for node_id, seq_len in seq_lens.items():
            node_label = ', '.join((node_id, str(seq_len)))
            node_size = calc_node_size(seq_len)
            fout.add_node(node_id,
                    {'label': node_label, 'fixedsize': 'true',
                        'height': str(node_size / 2), 'width': str(node_size)}
                    )

        # Write edges.
        for edge in edges:
            # write a edge in dot file.
            edge_colors = ['red;0.5' if ele < 0 else 'green;0.5' for ele in edge]
            edge_labels = [str(abs(ele)) for ele in edge]
            fout.add_edge(*edge_labels, {"color":':'.join(edge_colors)})


def print_help():
    body = '[-c] <input_file> [output_file]'
    print("python3 {} {}".format(__file__, body))

def main():
    interface = 'hc'
    options, args = getopt.getopt(sys.argv[1:], interface)
    to_compile_dot = False
    for option, value in options:
        if option == '-h':
            print_help()
            sys.exit()
        elif option == '-c':
            to_compile_dot = True
        else:
            print_help()
            sys.exit()

    input_file = args[0]
    if len(args) > 1:
        output_file = args[1]
    else:
        output_file = input_file + '.dot'

    write_dot_file(input_file, output_file)
    print('Generating dot file...')

    if to_compile_dot:
        print('Generating JPG file...')
        jpg_file = output_file.split('.', 1)[0] + '.jpg'
        subprocess.run(['dot', '-Tjpg', '-o', jpg_file, output_file])


if __name__ == '__main__':
    main()

