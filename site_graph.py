import difflib


class Site:

    def __init__(self, id_):
        self.id = id_
        self.children = []

    def add_child(self, child_site, interval, nodes_path, contamination_info=None):
        """
        Arguments: 
            child_site(Site): site that can be reach form current site.
            interval(int): Steps to take from the first base of start site to the
                first base of end site.    
            nodes_path(list[assembly_graph.Node]): length 1 or >1, when 1, it's it means
                there is no jump between nodes. When >2,  there are jump(s), the first and
                last node are nodes before and after jump(s)
            contamination_info(list[bool]): each element correspond to a ele in `node_path`, 
                indicate whether that node is contaminated, if it is, then that node's seq is replace by N.
        """
        if not contamination_info:
            contamination_info = [False] * len(nodes_path)
        self.children.append((child_site, interval, nodes_path, contamination_info))

    def __str__(self):
        return "{{Site {}}}".format(self.id)

    def __repr__(self):
        return str(self)

def simplify_site_graph(sites):
    num_edge_original = 0
    num_edge_later = 0
    for site in sites.values():
        num_edge_original += len(site.children)
        interval_endsite = {}
        new_children = []
        path_matchers = {}
        for child_site, interval, nodes_path, contamination in site.children:
            if (interval, child_site) not in path_matchers:
                path_matchers[(interval, child_site)] = (difflib.SequenceMatcher(None, nodes_path), nodes_path, contamination)
            else:
                matcher, template_nodes_path, contamination = path_matchers[(interval, child_site)]
                matcher.set_seq2(nodes_path)
                match_result = matcher.get_matching_blocks()
                index = 0
                for start_index, _, length in match_result:
                    while index < start_index:
                        contamination[index] = True
                        index += 1
                    index += length
        site.children = []
        for (interval, child_site), (_, nodes_path, contamination) in path_matchers.items():
            site.children.append((child_site, interval, nodes_path, contamination))
        num_edge_later += len(site.children)
    print('Simplified site graph: {} -> {}'.format(num_edge_original, num_edge_later))
    return sites

def read_file(file_name, nodes):
    fin = open(file_name)
    memory = {}
    sites = {}

    for line in filter(lambda x: len(x) > 3 and not x.startswith('#'), fin):
        tokens = line.split(' : ')
        start_site_id, end_site_id, interval = tokens[0].split()
        interval = int(interval)
        paths, is_contaminated_infos = [], []
        for node_id in tokens[1].rstrip().split():
            paths.append(nodes[node_id.rstrip('_')])
            is_contaminated_infos.append(node_id.endswith('_'))
        if start_site_id not in sites:
            sites[start_site_id] = Site(start_site_id)
        memo_info = (end_site_id, interval, paths, is_contaminated_infos)
        try:
            memory[start_site_id].append(memo_info)
        except KeyError:
            memory[start_site_id] = [memo_info]
    for start_site_id, children_info in memory.items():
        start_site = sites[start_site_id]
        for end_site_id, interval, paths, is_contaminated_infos in children_info:
            try:
                end_site = sites[end_site_id]
            except KeyError:
                end_site = Site(end_site_id)
                sites[end_site_id] = end_site
            start_site.addChild(end_site, interval, paths, is_contaminated_infos)

    fin.close()
    print('Load {} sites.'.format(len(sites)))
    return sites

def _test_read_file():
    input_file = 'assembly_graph.siteGraph'
    fastg_file_ = 'assembly_graph.fastg'
    import fastg_file
    nodes = fastg_file.build_assembly_graph(fastg_file_, 127)
    sites = read_file(input_file, nodes)
    for site in sites.values():
        print(site.id)
        print([(ele[0], ele[1]) for ele in site.children])

def write_file(file_name, sites, comment_list):
    fout = open(file_name, 'w')

    # Write commons
    for ele in comment_list:
        fout.write(''.join(('# ', ele, '\n')))

    sites_list = list(sites.values())
    sites_list.sort(key=lambda x: int(x.id.rstrip('r')))
    
    for i in range(0, len(sites_list), 2):
        if sites_list[i].id.endswith('r'):
            tmp = sites_list[i: i + 2]
            tmp.reverse()
            sites_list[i: i + 2] = tmp

    for site in sites_list:
        site.children.sort(key=lambda x: x[1])
        site.children.sort(key=lambda x: int(x[0].id.rstrip('r')))
        for child_site, interval, path, contamination in site.children:
            path_str = ' '.join((ele.uid + ('_' if contaminated else '') for (ele, contaminated) in zip(path, contamination)))
            fout.write('{} {} {} : {}\n'.format(site.id, child_site.id, interval, path_str))
        fout.write('\n')

    fout.close()


if __name__ == '__main__':
    _test_read_file()
