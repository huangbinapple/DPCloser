import sys
import getopt
import gap_info_file
import fastg_file
import sitegraph_builder
import site_graph
import spades_util


def print_help():
    body = '-f <fastg_file> -o <overlap size> [-s <site graph file>] <gap_info_file> <gap_path> <gap_seq>'
    print("python3 {} {}".format(__file__, body))

def get_site_by_index(index, site_position_index):
    positions = sorted(list(site_position_index.keys()))
    return site_position_index[positions[index]]

def get_node_id_from_long_name(long_name):
    short_name, _ = spades_util.read_long_name(long_name)
    id_, _, _, is_reversed = spades_util.read_short_name(short_name)
    return id_ + ('r' if is_reversed else '')

def main():
    # unmodified_reference_cmap = None
    fastg_file_name = None
    site_graph_file_name = None
    interface = 'hf:o:s:'
    options, args = getopt.getopt(sys.argv[1:], interface)
    for option, value in options:
        if option == '-h':
            print_help()
            sys.exit()
        elif option == '-f':
            fastg_file_name = value
        elif option == '-s':
            site_graph_file_name = value
        elif option == '-o':
            overlap_len = int(value)
    gap_info_file_name, gap_path_file, gap_seq_file = args

    # Read assembly graph.
    nodes = fastg_file.build_assembly_graph(
        fastg_file_name, overlap_len)

    if site_graph_file_name:
        # Reads site graph from file.
        sites = site_graph.read_file(site_graph_file_name, nodes)
        _, site_position_index = sitegraph_builder.build_site_graph(
            nodes, mode=1)
    else:
        # Build site graph form assembly graph, and write site graph.
        sites, site_position_index = \
            sitegraph_builder.build_site_graph(nodes)
        sites = site_graph.simplify_site_graph(sites)
        site_graph_file_name = fastg_file_name.rsplit('.')[0] +\
             '.sitegraph'
        site_graph.write_file(site_graph_file_name, sites)

    # Find start site and end site.
    gaps = gap_info_file.read_file(gap_info_file_name)

    for gap in gaps:
        start_node_id = get_node_id_from_long_name(gap.start_node_id)
        end_node_id = get_node_id_from_long_name(gap.end_node_id)
        start_site = get_site_by_index(gap.start_site_index, 
            site_position_index[nodes[start_node_id]])
        end_site = get_site_by_index(gap.end_site_index, 
            site_position_index[nodes[end_node_id]])
        print(start_site.id, end_site.id)


if __name__ == '__main__':
    main()
        