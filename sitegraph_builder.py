import sys
import getopt
import time

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

def infer_site_position_in_reverse_complement(position, seq_length, key_word_length, shift=0):
    return seq_length - key_word_length - position - shift

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

def add_sites_to_node_pair(node, site_position_indexs, min_id_avalable, nodes, sites, shift=0):
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
            site_position_ic = infer_site_position_in_reverse_complement(site_position, node.length, len(BIONANO_SITE[0]), shift)
            
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

ALLOWED_REPEAT_NUM = 2

def build_site_graph(nodes, shift=0, mode=0, max_interval_len=100000):
    print('len of nodes:', len(nodes))
    sites = {}
    site_position_index = {}
    node_id = 0
    for node in filter(lambda x: x.uid[-1] != 'r', sorted(nodes.values(), key=lambda x: x.uid)):
        node_ic = get_node_ic(node, nodes)
        num_sites_added = add_sites_to_node_pair(node, site_position_index, node_id, nodes, sites, shift)
        node_id += num_sites_added
        # if num_sites_added:
            # print('add {} site(s) in node {}, {}'.format(num_sites_added, node.uid, node_ic.uid))

    if mode == 1:
        # only return index:
        return None, site_position_index
    # Add site graph edges that jump between nodes.
    print('Number of sites: {}'.format(len(sites)))
    print('Number of nodes: {}'.format(len(nodes)))
    print('Numner of nodes contain site: {}'.format(len(site_position_index)),)
    print('Great percentage:', round(len(site_position_index) / len(nodes) * 100, 2))
    for node, site_positions in site_position_index.items():
        # if node.uid != '252560r':
            # continue
        max_position = max(site_positions.keys())
        max_position_site = site_positions[max_position]
        print('Add child to: {}'.format(max_position_site))
        interval = 0
        node_path = []
        sub_intervals = []
        # print('node.length:', node.length)
        # print('max_position:', max_position)
        sub_interval = node.length - max_position - 1  # Steps to take from `max_position` to end of node.
        frontier = [(node, sub_interval)]
        # debug = False
        while frontier:
            ele = frontier.pop()
            if type(ele) == tuple:
                frontier.append(None)
                current_node, sub_interval = ele
                # print('sub_interval:', sub_interval)
                interval += sub_interval
                node_path.append(current_node)
                sub_intervals.append(sub_interval)
                # if current_node.uid == '252560r':
                    # debug = True
                    # print('haha')

                for child_node, overlap_len in current_node.children:
                    if node_path.count(child_node) >= ALLOWED_REPEAT_NUM or interval > max_interval_len:
                        continue
                    try:
                        to_continue = False
                        child_site_positions_index = site_position_index[child_node]
                        child_site_positions = list(child_site_positions_index.items())
                        child_site_positions.sort()
                        min_child_site_position = overlap_len - len(BIONANO_SITE[0])  # Child_site's position should > this number.
                        for child_site_position, child_site in child_site_positions:
                            if child_site_position > min_child_site_position:
                                # if debug:
                                    # print(max_position_site.id, child_site.id, [ele.uid for ele in node_path] + [child_node.uid])
                                # print('interval:', interval)
                                # print('overlap_len:', overlap_len)
                                # print('last_site_positoin:', child_site_position)
                                # print('last sub_interval:', interval + child_site_position - overlap_len + 1)
                                max_position_site.add_child(child_site, interval + child_site_position - overlap_len + 1,
                                        node_path + [child_node])
                                # print('Add a child to {}'.format(max_position_site.id), '*' * (len(node_path) + 1))
                                to_continue = True
                                break
                        if to_continue:
                            continue
                        raise KeyError('This is a Fake KeyError, just mean we should continue searchiing for child site.')
                    except KeyError:
                        frontier.append((child_node, child_node.length - overlap_len))
            elif ele is None:
                last_node, sub_interval_ = node_path.pop(), sub_intervals.pop()
                # if last_node.uid == '253562':
                    # debug = False
                    # print('Oh')
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
        # print(site, positions, 'C:', child_site_id_intervals_)
        if site.id in child_site_ids:
            num_site_contain_self_as_child += 1
            print('!!!')
        count[len(positions)] += 1


def printHelpMessage():
    body = '[-h] [-s]  <-i fastg file> <-m max_interval_len> <-l overlap_len> <-o site_graph file>'
    print('python3 {} {}'.format(__file__, body))


def main():
    input_file = ''
    output_file = ''
    overlap_len = None
    shift_len = 0
    to_simplify = False
    max_interval_len = None
    options, args = getopt.getopt(sys.argv[1:], 'm:k:i:l:o:hs')
    for option, value in options:
        if option == '-i':
            input_file = value
        elif option == '-o':
            output_file = value
        elif option == '-l':
            overlap_len = int(value)
        elif option == '-k':
            shift = int(value)
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
    print('overlap_len:', overlap_len)
    nodes = fastg_file.build_assembly_graph(input_file, overlap=overlap_len)
    tick = time.time()
    sites, site_position_index = build_site_graph(nodes, shift=shift)
    tock = time.time()

    if to_simplify:
        print('Simplifying site graph...')
        sites = site_graph.simplify_site_graph(sites)
    # print('{} sites created on {} nodes'.format(len(sites), len(site_position_index)))
    tock2 = time.time()

    # debug_site = sites['648r']
    # print('site', debug_site.id, 'has {} children.'.format(len(debug_site.children)))
    # for child, interval, nodes_path, _ in debug_site.children:
        # print(child.id, interval, [ele.uid for ele in nodes_path])

    site_graph.write_file(output_file, sites,
            [' '.join(sys.argv),
             'Number of sites: {}'.format(len(sites)),
             'Number of nodes: {}'.format(len(nodes)),
             'Numner of nodes contain site: {}'.format(len(site_position_index)),
             'Max interval len: {}'.format(max_interval_len),
             'Allowed repeat num: {}'.format(ALLOWED_REPEAT_NUM),
             'Time used to build graph: {} seconds'.format(round(tock - tick)),
             'Time used to simplify graph: {} seconds'.format(round(tock2 - tock))]
            )

if __name__ == '__main__':
    main()
