"""This file is used to locate sites' position on spades EDGE(assembly graph node) after they have been
located on spades NODE(contig)"""

import sys
import getopt

import paths_file
import gap_info_file
import fastg_file 
import sitegraph_builder


def get_site_position(nodes, node_id_list, site_position_indexs, site_rank, len_overlap):
    sites_in_order = []
    for i, node_id in enumerate(node_id_list):
        node = nodes[node_id]
        try:
            site_positions = site_position_indexs[node].items()
        except KeyError:
            continue
        if i != 0:
            site_positions = filter(lambda x: x[0] >= len_overlap, site_positions)
        site_positions = sorted(list(site_positions))
        # print('haha')
        num_site = len(site_positions)
        # print('num_site:', num_site)
        if not site_positions:
            continue
        sites_in_order_sub = list(zip(*site_positions))[1]
        # print(list(sites_in_order_sub))
        sites_in_order.extend(list(zip([node] * num_site, sites_in_order_sub)))
    # for i, ele in enumerate(sites_in_order):
        # print(i, ele)
    result_node, result_site = sites_in_order[site_rank]
    # print(result_node, result_site)
    result_site_index, result_site_position = None, None
    # print(site_position_indexs[result_node])
    for i, (site_position, site_) in enumerate(sorted(site_position_indexs[result_node].items())):
        if site_ == result_site:
            result_site_index, result_site_position = i, site_position
            break
    return result_node.uid, result_site_index, result_site_position


def transform_position(original_node_id, original_site_index, nodes, paths, site_position_indexs, len_overlap):
    node_id_list = paths[original_node_id]
    result_node_id, result_site_index, result_site_position = get_site_position(
        nodes, node_id_list, site_position_indexs, original_site_index, len_overlap=OVERLAP)
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
        gap.start_node_id, gap.start_site_index, gap.start_site_position =\
            transform_position(gap.start_node_id, gap.start_site_index, nodes, paths, site_position_indexs, OVERLAP)
        gap.end_node_id, gap.end_site_index, gap.end_site_position =\
            transform_position(gap.end_node_id, gap.end_site_index, nodes, paths, site_position_indexs, OVERLAP)
    # Write gaps.
    cmd = ' '.join(sys.argv)
    gap_info_file.write_file(gaps, output_file_name, comments=[cmd])


if __name__ == '__main__':
    main()

