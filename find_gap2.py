"""This file is used to locate sites' position on spades EDGE(assembly graph node) after they have been
located on spades NODE(contig)"""

import sys
import getopt
import bisect

import paths_file
import gap_info_file
import fastg_file 
import sitegraph_builder


def get_site_position(nodes, node_paths, site_position_indexs, site_position, len_overlap, debug=0):
    # print('node_paths:', node_paths)
    # print('site_position:', site_position)
    current_node_index = 0
    current_node = nodes[node_paths[current_node_index]]
    while site_position > current_node.length:
        next_node_id = node_paths[current_node_index + 1]
        # print('next_node_id:', next_node_id)
        child_index = 0
        while current_node.children[child_index][0].uid != next_node_id:
            child_index += 1
        next_node, overlap = current_node.children[child_index]
        site_position -= (current_node.length - overlap)
        current_node = next_node
        current_node_index += 1
    site_position -= 1  # Convert 1-index to 0-index.
    # print(current_node, site_position)

    site_positions = site_position_indexs[current_node]
    if debug:
        print('current_node_uid: ', current_node.uid)
        for k, v in sorted(list(site_positions.items())):
            print(k, v)
    try:
        result_site = site_positions[site_position]
        site_index = sorted(list(site_positions.keys())).index(site_position)

    except KeyError:
        positions_in_order = sorted(list(site_positions.keys()))
        upper_index = bisect.bisect(positions_in_order, site_position)
        lower_index = upper_index - 1
        # print('positions_in_order:', positions_in_order)
        # print('site_info:', site_info)
        # print('site_position:', site_position)
        # print('site_position_in_contig:', site_position_in_contig)
        # print('current_node:', current_node)
        # print('neighbors:', positions_in_order[lower_index:upper_index+1])
        # assert abs(sum(positions_in_order[lower_index:upper_index+1]) / 2 - site_position) < 1
        # adjusted_site_position = positions_in_order[lower_index if down_when_confused else upper_index]
        site_position = positions_in_order[upper_index]
        result_site = site_positions[site_position]
        site_index = upper_index

    return current_node.uid, site_index, site_position

def transform_position(original_node_id, original_site_position, nodes, paths, site_position_indexs, len_overlap, debug=0):
    node_id_list = paths[original_node_id]
    result_node_id, result_site_index, result_site_position = get_site_position(
        nodes, node_id_list, site_position_indexs, original_site_position, len_overlap=OVERLAP, debug=debug)
    return result_node_id, result_site_index, result_site_position

def print_help():
    body = '-p <path file> -g <fastg file> <gap info file> <output file>'
    print("python3 {} {}".format(__file__, body))

OVERLAP = 77
def main():

    paths_file_name = ''
    gap_info_file_name = ''
    fastg_file_name = ''
    interface = 'hp:g:'
    options, args = getopt.getopt(sys.argv[1:], interface)
    for option, value in options:
        if option == '-h':
            print_help()
            sys.exit()
        elif option == '-p':
            paths_file_name = value
        elif option == '-g':
            fastg_file_name = value
    gap_info_file_name, output_file_name = args

    # paths_file_name = '/home/huangbin/simulation_lab/Bacteria/E.coli_PL/spades.dir/contigs.paths'
    # gap_info_file_name = '/home/huangbin/simulation_lab/Bacteria/E.coli_PL/hybrid_scaffold_keep_11.dir/hybrid_scaffolds/tmp1.txt'
    # fastg_file_name = 'datasets/assembly_graph.fastg'
    paths = paths_file.read_file(paths_file_name)
    gaps = gap_info_file.read_file(gap_info_file_name)
    nodes = fastg_file.build_assembly_graph(fastg_file_name, overlap=OVERLAP)
    _, site_position_indexs = sitegraph_builder.build_site_graph(nodes, mode=1)
    # Transform gaps.
    for gap in gaps:
        debug = 0
        if gap.start_node_id.startswith('NODE_2_'):
            debug = 1
        if debug == 1:
            print('DEBUG:', gap.start_node_id, gap.start_site_index)
        gap.start_node_id, gap.start_site_index, gap.start_site_position =\
            transform_position(gap.start_node_id, gap.start_site_position, nodes, paths, site_position_indexs, OVERLAP, debug=debug)
        gap.end_node_id, gap.end_site_index, gap.end_site_position =\
            transform_position(gap.end_node_id, gap.end_site_position, nodes, paths, site_position_indexs, OVERLAP, debug=debug)
        if debug == 1:
            print('DEBUG:', gap.start_node_id, gap.start_site_index)
    # Write gaps.
    cmd = ' '.join(sys.argv)
    gap_info_file.write_file(gaps, output_file_name, comments=[cmd])


if __name__ == '__main__':
    main()

