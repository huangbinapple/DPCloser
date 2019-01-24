import sys
import getopt

import fastg_file
import itertools
from site_graph import Site
import site_graph


BIONANO_SITE = 'GCTCTTC', 'GAAGAGC'

def _find_all_helper(string, keyword, result):
    hitIndex = string.find(keyword)
    while hitIndex >= 0:
        result.append(hitIndex)
        hitIndex = string.find(keyword, hitIndex+1)

def find_all(string, keywords):
    result = []
    for keyword in keywords:
        _find_all_helper(string, keyword, result)
    result.sort()
    return result

def test_find_all():
    # BIONANO_SITE = 'GA', 'TC'
    # string = "AATCCCCCCCTCACCGCCATATTTAAAGAAGAGCTCGTACGAAAGTACGGGCTTTTTTTTCGTATATTGCACACACCGGGGGGATGAGAAGCCCCGACCGGGGTTCGACAACTGGCGACAGCCAGTTGGACAGACCGTGAACGCAGTGAGCGGGCTGCCCGCAGGGCGAGCGAAGCGAGTCAATCCCCCCCTCACCGCCATATTTAAAGAAGAGCTCGTACGAAAGTACGGGCTTTTTTTTCGTATATTGCACACACCGGG"
    string_ic = "CCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGCTCTTCTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCCCGCTCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGAACCCCGGTCGGGGCTTCTCATCCCCCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGCTCTTCTTTAAATATGGCGGTGAGGGGGGGATT"
    print(find_all(string_ic, BIONANO_SITE))

def infer_site_position_in_reverse_complement(position, seq_length, key_word_length):
    return seq_length - key_word_length - position

def infer_site_position_in_child(position, parent_length, overlap_length):
    return position - (parent_length - overlap_length)

def is_site_in_overlap(position, node_len, over_lap_len):
    return position >= (node_len - over_lap_len)

def get_node_ic(node, nodes):
    """Get the inverse complement node of `node` in `nodes`."""
    result_node_id = node.uid[:-1] if node.uid.endswith('r') else node.uid + 'r'
    return nodes[result_node_id]

def get_site_ic(site, sites):
    """Get the inverse complement node of `site` in `sites`."""
    result_site_id = site.id[:-1] if site.id.endswith('r') else site.id + 'r'
    return sites[result_site_id]

def add_sites_to_node_pair(node, site_position_indexs, min_id_avalable, nodes, sites):
    node_ic = get_node_ic(node, nodes)

    # Find positions, if there is none, return.
    site_positions = find_all(node.seq, BIONANO_SITE)
    if site_positions:
        position_index, position_index_ic = {}, {}
        site_position_indexs[node] = position_index
        site_position_indexs[node_ic] = position_index_ic

        for site_position in site_positions:
            new_site_id = str(min_id_avalable)
            new_site = Site(new_site_id)
            new_site_ic = Site(new_site_id + 'r')
            site_position_ic = infer_site_position_in_reverse_complement(site_position, node.length, len(BIONANO_SITE[0]))
            
            # Register new site(_ic) in `position_index`(_ic) & `sites`.
            position_index[site_position] = new_site
            position_index_ic[site_position_ic] = new_site_ic
            sites[new_site_id] = new_site
            sites[new_site_id + 'r'] = new_site_ic

            min_id_avalable += 1

        # Add site graph edges within one node.
        for node_, position_index_ in ((node, position_index), (node_ic, position_index_ic)):
            positions, sites_in_order = zip(*sorted(list(position_index_.items())))
            for i in range(len(sites_in_order) - 1):
                sites_in_order[i].add_child(sites_in_order[i + 1], positions[i + 1] - positions[i], [node_])
    return len(site_positions)

def build_site_graph(nodes, mode=0, max_interval_len=100000):
    print('len of nodes:', len(nodes))
    sites = {}
    site_position_index = {}
    node_id = 0
    for node in filter(lambda x: x.uid[-1] != 'r', sorted(nodes.values(), key=lambda x: x.uid)):
        node_ic = get_node_ic(node, nodes)
        num_sites_added = add_sites_to_node_pair(node, site_position_index, node_id, nodes, sites)
        node_id += num_sites_added
        if num_sites_added:
            print('add {} site(s) in node {}, {}'.format(num_sites_added, node.uid, node_ic.uid))

    if mode == 1:
        # only return index:
        return None, site_position_index
    # Add site graph edges that jump between nodes.
    for node, site_positions in site_position_index.items():
        max_position = max(site_positions.keys())
        max_position_site = site_positions[max_position]
        interval = 0
        allowed_repeat_num = 3
        node_path = []
        sub_intervals = []
        sub_interval = node.length - max_position - 1  # Steps to take from `max_position` to end of node.
        frontier = [(node, sub_interval)]
        while frontier:
            ele = frontier.pop()
            if type(ele) == tuple:
                frontier.append(None)
                current_node, sub_interval = ele
                interval += sub_interval
                node_path.append(current_node)
                sub_intervals.append(sub_interval)
                for child_node, overlap_len in current_node.children:
                    if node_path.count(child_node) >= allowed_repeat_num or interval > max_interval_len:
                        continue
                    try:
                        to_continue = False
                        child_site_positions_index = site_position_index[child_node]
                        child_site_positions = list(child_site_positions_index.items())
                        child_site_positions.sort()
                        min_child_site_position = overlap_len - len(BIONANO_SITE[0])  # Child_site's position should > this number.
                        for child_site_position, child_site in child_site_positions:
                            if child_site_position > min_child_site_position:
                                max_position_site.add_child(child_site, interval + child_site_position - overlap_len + 1,
                                        node_path + [child_node])
                                print('Add a child to {}'.format(max_position_site.id), '*' * (len(node_path) + 1))
                                to_continue = True
                                break
                        if to_continue:
                            continue
                        raise KeyError('This is a Fake KeyError, just mean we should continue searchiing for child site.')
                    except KeyError:
                        frontier.append((child_node, child_node.length - overlap_len))
            elif ele is None:
                last_node, sub_interval_ = node_path.pop(), sub_intervals.pop()
                interval -= sub_interval_
            
    return sites, site_position_index

def _test_build_site_graph():
    input_file = 'assembly_graph.fastg'
    # input_file = 'test.fastg'
    overlap = 127
    nodes = fastg_file.build_assembly_graph(input_file, overlap)
    # nodes = fastg_file.build_assembly_graph(input_file, overlap=2)
    sites, site_position_index = build_site_graph(nodes)
    num_position_in_index = sum((len(ele) for ele in site_position_index.values()))
    site_ids = [int(site.id) for site in sites.values() if not site.id.endswith('r')]
    site_ids.sort()

    site_positions = list(site_position_index.items())
    site_positions.sort(key=lambda x: int(x[0].uid.rstrip('r')))

    site_index = {}
    count = [0] * 10
    num_site_contain_self_as_child = 0
    for node, positions in site_position_index.items():
        for position, site in positions.items():
            if site not in site_index:
                site_index[site] = [(node.uid, position)]
            else:
                site_index[site].append((node.uid, position))
    for v in site_index.values():
        v.sort(key=lambda x: int(x[0].rstrip('r')))
    site_index_items = list(site_index.items())
    site_index_items.sort(key=lambda x: x[1][0])
    for site, positions in site_index_items:
        child_sites, intervals, node_paths = list(zip(*site.children)) if site.children else ([], [], [])
        child_site_ids = [ele.id for ele in child_sites]
        node_ids = [[ele.uid for ele in node_path] for node_path in node_paths]
        child_site_id_intervals = list(zip(child_site_ids, intervals))
        child_site_id_intervals_ = sorted(child_site_id_intervals)
        print(site, positions, 'C:', child_site_id_intervals_)
        if site.id in child_site_ids:
            num_site_contain_self_as_child += 1
            print('!!!')
        count[len(positions)] += 1

    print('count:', count)
    print('len of sites:', len(sites))
    print('len of site position index:', len(site_position_index))
    print('number of position:', num_position_in_index)
    print('number of self loop:', num_site_contain_self_as_child)


def printHelpMessage():
    body = '[-h] <-i fastg file> <-m max_interval_len> <-l overlap_len> <-o site_graph file>'
    print('python3 {} {}'.format(__file__, body))


def main():
    input_file = ''
    output_file = ''
    overlap_len = None
    to_simplify = False
    max_interval_len = None
    options, args = getopt.getopt(sys.argv[1:], 'm:i:l:o:hs')
    for option, value in options:
        if option == '-i':
            input_file = value
        elif option == '-o':
            output_file = value
        elif option == '-l':
            overlap_len = int(value)
        elif option == '-m':
            max_interval_len = int(value)
        elif option == '-s':
            to_simplify = True
        elif option == '-h':
            printHelpMessage()
            sys.exit()
        else:
            printHelpMessage()
            sys.exit()

    nodes = fastg_file.build_assembly_graph(input_file, overlap=overlap_len)
    sites, site_position_index = build_site_graph(nodes)
    if to_simplify:
        print('Simplifying site graph...')
        sites = site_graph.simplify_site_graph(sites)
    print('{} sites created on {} nodes'.format(len(sites), len(site_position_index)))

    site_graph.write_file(output_file, sites,
            [' '.join(sys.argv),
             'Number of sites: {}'.format(len(sites)),
             'Number of nodes: {}'.format(len(nodes)),
             'Numner of nodes contain site: {}'.format(len(site_position_index))]
            )

if __name__ == '__main__':
    main()
