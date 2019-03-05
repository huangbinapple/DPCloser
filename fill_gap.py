import sys
import getopt
import pickle
import subprocess
from multiprocessing import Pool
import gap_info_file
import fastg_file
import sitegraph_builder
import site_graph
import spades_util
import find_path_dp
sys.setrecursionlimit(10000)


def get_site_by_index(index, site_position_index):
    positions = sorted(list(site_position_index.keys()))
    position = positions[index]
    return position, site_position_index[position]

def get_node_id_from_long_name(long_name):
    short_name, _ = spades_util.read_long_name(long_name)
    id_, _, _, is_reversed = spades_util.read_short_name(short_name)
    return id_ + ('r' if is_reversed else '')

def construct_cmd(start_site_id, end_site_id, start_site_position,
    end_site_position, graph_pickle, intervals, gap_seq_id):
    exe_list = []
    if start_site_id:
        exe_list.extend(['-s', start_site_id])
    if end_site_id:
        exe_list.extend(['-e', end_site_id])
    if start_site_position:
        exe_list.extend(['-a', str(start_site_position)])
    if end_site_position:
        exe_list.extend(['-b', str(end_site_position)])
    exe_list.extend(['-n', gap_seq_id])
    exe_list.extend(['-g', graph_pickle_file])
    exe_list.extend(['-i', ','.join(map(str, intervals))])
    exe_list.extend(['-o', seq_fa_file])
    return exe_list

def find_path_dp_process(exe_tuple):
    mol_name, cmd_list = exe_tuple
    with open(log_dir + '/' + mol_name + '.log', 'w') as stdout:
        print('processing {}...'.format(mol_name))
        res = subprocess.run(['python3', find_path_script] + cmd_list,
            stdout=stdout)
        print(mol_name, 'SUCCESS' if res.returncode == 0 else 'FAIL')


def print_help():
    body = '-f <fastg_file> -o <overlap size> -x <find path script dir> [-s <site graph file>] <gap info file> [work_dir] [task_name]'
    '<gap_info_file> <graph pickle file> <log_dir> <gap_fa>'
    print("python3 {} {}".format(__file__, body))

def main():
    # unmodified_reference_cmap = None
    global find_path_script
    fastg_file_name = None
    site_graph_file_name = None
    find_path_script = 'find_path_dp.py'
    work_dir = None
    task_name = None
    n_thread = 1
    is_node_id_processed = 0
    interface = 'hf:o:n:x:s:m:'
    options, args = getopt.getopt(sys.argv[1:], interface)
    for option, value in options:
        if option == '-h':
            print_help()
            sys.exit()
        elif option == '-f':
            fastg_file_name = value
        elif option == '-s':
            site_graph_file_name = value
        elif option == '-x':
            find_path_script = value
        elif option == '-n':
            n_thread = int(value)
        elif option == '-o':
            overlap_len = int(value)
        elif option == '-m':
            is_node_id_processed = int(value)
    
    gap_info_file_name = args[0]
    if len(args) == 1:
        pass
    if len(args) == 2:
        task_name = args[1]
    elif len(args) == 3:
        work_dir, task_name = args[1:3]
    else:
        print_help()
        sys.exit()
    if task_name is None:
        task_name = fastg_file_name.rsplit('.', 1)
    global log_dir, seq_fa_file, graph_pickle_file
    log_dir = work_dir + '/' + task_name + '_log.dir'
    seq_fa_file = work_dir + '/' + task_name + '.fa'
    graph_pickle_file = work_dir + '/' + task_name + '_graph.pickle'
    subprocess.run(('rm', '-rf', log_dir))
    subprocess.run(('mkdir', log_dir))

    # Read assembly graph.
    nodes = fastg_file.build_assembly_graph(
        fastg_file_name, overlap_len)

    if site_graph_file_name:
        # Reads site graph from file.
        sites = site_graph.read_file(site_graph_file_name, nodes)
        _, site_position_index = sitegraph_builder.build_site_graph(
            nodes, mode=1)
        for position_index in site_position_index.values():
            for position in position_index:
                original_site = position_index[position]
                position_index[position] = sites[original_site.id]
    else:
        # Build site graph form assembly graph, and write site graph.
        sites, site_position_index = \
            sitegraph_builder.build_site_graph(nodes)
        sites = site_graph.simplify_site_graph(sites)
        site_graph_file_name = fastg_file_name.rsplit('.')[0] +\
             '.sitegraph'
        site_graph.write_file(site_graph_file_name, sites)

    # Write nodes and sites to pickle file.
    with open(graph_pickle_file, 'wb') as fout:
        pickle.dump(nodes, fout, -1)
        pickle.dump(sites, fout, -1)

    # Find start site and end site.
    gaps = gap_info_file.read_file(gap_info_file_name)

    exe_tuples = []
    for gap in gaps:
        if is_node_id_processed:
            start_node_id = gap.start_node_id
            end_node_id = gap.end_node_id
        else:
            start_node_id = get_node_id_from_long_name(gap.start_node_id)
            end_node_id = get_node_id_from_long_name(gap.end_node_id)
        start_site_position, start_site = get_site_by_index(
            gap.start_site_index, 
            site_position_index[nodes[start_node_id]])
        end_site_position, end_site = get_site_by_index(
            gap.end_site_index, 
            site_position_index[nodes[end_node_id]])
        print(start_site.id, end_site.id)

        gap_seq_id = '-'.join((start_node_id, end_node_id))
        exe_tuples.append((gap_seq_id,
            construct_cmd(start_site.id, end_site.id,
                start_site_position, end_site_position,
                graph_pickle_file, gap.intervals, gap_seq_id)))

    with Pool(n_thread) as p:
        p.map(find_path_dp_process, exe_tuples)



if __name__ == '__main__':
    main()
                
