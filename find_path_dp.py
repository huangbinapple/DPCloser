import math
import pickle
import sys
import getopt
import statistics
from copy import copy
import numpy as np
import site_graph
import time
import difflib


MAX_CONTINUE_FP = 1
MAX_CONTINUE_FN = 1.9
print('MAX_CONTINUE_FP:', MAX_CONTINUE_FP)
print('MAX_CONTINUE_FN:', MAX_CONTINUE_FN)

def get_latest_row(matrix, end_index):
    """
    Argument:
        matrix (numpy.array, size: m x any_col): matrix we should retrieve data from.
        end_index (int): smalleset index of row which was NOT retrived.
    Return:
        (numpy.array, size: n x any_col. n <= MAX_CONTINUE_FP + 1 + 1)
    """
    start_index = max(0, end_index - MAX_CONTINUE_FP - 1)
    return matrix[start_index:end_index, :]

def get_latest_interval_sums(intervals, end_index):
    """
    Return:
        (numpy.array, size: n x 1, n <= MAX_CONTINUE_FP + 1)
    """
    start_index = max(0, end_index - MAX_CONTINUE_FP - 1)
    intervals_ = intervals[start_index:end_index]
    sum_ = 0
    result = []
    for interval in reversed(intervals_):
        sum_ += interval
        result.append(sum_)
    return np.array(list(reversed(result))).reshape(-1, 1)

FACTORIAL_N = [1]
for i in range(1, MAX_CONTINUE_FP + 1):
    FACTORIAL_N.append(FACTORIAL_N[-1] * i)
FACTORIAL_N_INVERSE = [1 / ele for ele in FACTORIAL_N]

def factorial_inverse_helper(x):
    return FACTORIAL_N_INVERSE[x]
factorial_inverse = np.vectorize(factorial_inverse_helper)

"""
We assume in interval L, the #insertion X obey possion distribution, expectation: LAMBDA = L / AVERAGE_INTERVAL.
which is: (L / AVERAGE_INTERVAL) ^ x / x! * exp(-L / AVERAGE_INTERVAL)
"""
# AVERAGE_INTERVAL = 74200
# AVERAGE_INTERVAL = 3053000
AVERAGE_INTERVAL = 1526000
print('AVERAGE_INTERVAL:', AVERAGE_INTERVAL)
def insert_factor(interval, num_insert):
    """
    Argument:
        interval (numpy.array, size: n x 1, dtype: int): intervals on the bionano molecules.
        num_insert (numpy.array, size: n x 1, dtype: int): number of site insert happened in the intervas.
    Return:
        (numpy.array, size: n x 1, dtype: float64): score factor represent the probability of insertion.
    """
    lambda_ = interval / AVERAGE_INTERVAL
    return lambda_ ** num_insert * np.exp(-lambda_) * factorial_inverse(num_insert)

"""
We assume that when reference's length is L, the measure is a laplace distribution, expetaion: MU, SD: SIGMA.
which is: 1 / DECAY * exp(-abs(x - (MU + L)) / DECAY),
after emit common factor, it is: exp(-abs(x - (MU + L)) / DECAY).
DECAY = sqrt(SIGMA / 2)
"""
MU = 0
# MU = 293  # alter
SIGMA = 200
FIVE_SIGMA = 5 * SIGMA
# SIGMA = 846  # alter
DECAY_INVERSE = math.sqrt(2) / SIGMA
print('MU:', MU)
print('SIGMA:', SIGMA)
print('FIVE_SIGMA:', FIVE_SIGMA)
def similar_factor(length_reference, length_bionano):
    """
    Argument:
        length_reference (int): path length in the site graph.
        length_bionano (numpy.array, size: n x 1, dtype: int): intervals on bionano molecules. 
    Return:
        (numpy.array, size: n x 1, dtype: float64): score factor represent the similarity between measured and reference length.
    """
    error = np.abs(length_bionano - (MU + length_reference))
    error = np.minimum(error, FIVE_SIGMA)
    result = np.exp(- error * DECAY_INVERSE)
    return result

def _test_similar_factor():
    print('DECAY_INVERSE:', DECAY_INVERSE)
    result = []
    for length_ref, length_bionano in ((1133, np.array([25191])),):
        score = similar_factor(length_ref, length_bionano)
        print('{} ~ {} -> {}'.format(length_bionano, length_ref, score))
        result.append(score)
    # print('ratio:', result[0] / result[1])

MERGE_LEN = 4000
FN_RATE = .228
print('MERGE_LEN:', MERGE_LEN)
print('FN_RATE:', FN_RATE)
def prob_skip(interval):
    """
    The probability of a site is FN (so it will be skipped in propagate process), given its interval to formmer site.
    Argument:
        interval (int): interval to formmer site.
    Return:
        (float): Probability.
    """
    # return .1  #  For test
    assert interval >= 0
    result = FN_RATE
    if interval < MERGE_LEN:
        return (1 - FN_RATE) * (interval / MERGE_LEN - 1) ** 2 + FN_RATE
    else:
        return FN_RATE

def _test_prob_skip():
    for interval in (10000, 4897, 457, 355, 2411, 1523, 1647, 240, 4742):
        print(interval, prob_skip(interval))

def modify(prob_to_modify, tracker_to_modify, fingerprints_to_modify, prob_candidate, fingerprints_candidate, start_site, child_indexs):
    """
    Update prob and tracker matrix row.
    Arguments:
        prob_to_modify (numpy.array, size: n): probs before update.
        tracker_to_modify (numpy.array, size: n): trackers before update.
        fingerprints_to_modify (numpy.array, size: n): fingerprints before update.
        prob_candidate (numpy.array, size: m x n): probs possible to update `prob_to_update`.
        fingerprints_candidate (numpy.array, size: m x n): finger print of correspond probs in `prob_candidate`.
        start_site (site_graph.Site): A site whose table is used to modify other sites' table.
        child_indexs (list[int]): Route taken from `start_site` to the site whose table is modified.
    Return:
        None
    """
    probs = np.vstack((prob_candidate, prob_to_modify))
    fingerprints = np.vstack((fingerprints_candidate, fingerprints_to_modify))
    size_row, n_rank = probs.shape

    probs_flat = probs.flat
    tmp = sorted(zip(probs_flat, range(len(probs_flat)), fingerprints.flat), key=lambda x: x[0], reverse=True)

    prob_modified, tracker_modified, fingerprint_modified = np.zeros(n_rank), np.empty(n_rank, dtype='object'), np.empty(n_rank, dtype='uint64')
    recorded_fingerprints = set()
    counter = 0
    for prob, index, fingerprint in tmp:
        if fingerprint in recorded_fingerprints:
            continue
        prob_modified[counter] = prob
        fingerprint_modified[counter] = fingerprint
        index_row, index_col = index // n_rank, index % n_rank
        tracker_modified[counter] = tracker_to_modify[index_col] if index_row == size_row - 1 else\
            (start_site, (index_row - size_row + 1, index_col), copy(child_indexs))
        counter += 1
        recorded_fingerprints.add(fingerprint)
        if counter == n_rank:
            break
    prob_to_modify[:], tracker_to_modify[:], fingerprints_to_modify[:] = prob_modified, tracker_modified, fingerprint_modified

HASH_A = 4343534
HASH_B = 232423
def update_fingerprint(original_fingerprint, child_index):
    return (original_fingerprint + HASH_A) * (child_index + HASH_B)

CRITICAL_FN_FACTOR = FN_RATE ** (MAX_CONTINUE_FN)
print('CRITICAL_FN_FACTOR:', CRITICAL_FN_FACTOR)
def propagate(i, start_site, intervals, prob_matrixs, tracker_matrixs, fingerprint_matrixs, min_prob):
    """
    Alter tables of sites who's in the downstream of `start_site`.
    Arguments:
        i (int): A pointer on the bionano molecule.
        start_site (site_graphe.Site): Site whose table is used to alter downstream sites' table.
        intervals (list[float], length: m - 1): Intervals between adjacent sites on bionano molecule.
        prob_matrixs (dict[site_graph.Site: numpy.array of size m x n): All sites' prob table.
        tracker_matrixs (dict[site_graph.Site: numpy.array of size m x n): All sites' tracker table.
        fingerprint_matrixs (dict[site_graph.Site: numpy.array of size m x n): All sites' fingerprint table.
    Return:
        (Set): All sites reached from `start_site`.
    """
    result = set()
    prob_matrix, tracker_matrix, fingerprint_matrix = prob_matrixs[start_site], tracker_matrixs[start_site], fingerprint_matrixs[start_site]
    content = get_latest_row(prob_matrix, i)
    if content[:, 0].max() < min_prob:
        return set()
    fingerprints_ = get_latest_row(fingerprint_matrix, i)
    intervals_ = get_latest_interval_sums(intervals, i)
    fp_factor = insert_factor(intervals_, np.array(list(reversed(range(intervals_.shape[0])))).reshape(-1, 1))
    # content *= fp_factor  # This is a BUG
    content = content * fp_factor

    site_path, child_indexs, fn_factors, path_interval_sums, fingerprints_stack = [], [], [], [], []
    frontier = [(None, None, start_site)]
    while frontier:
        ele = frontier.pop()
        if type(ele) is tuple:
            child_index, current_interval, current_site = ele 
            frontier.append(None)

            site_path.append(current_site)
            child_indexs.append(child_index)
            fn_factors.append(fn_factors[-1] * prob_skip(current_interval) if fn_factors else 1)
            path_interval_sums.append(path_interval_sums[-1] + current_interval if path_interval_sums else 0)
            fingerprints_stack.append(update_fingerprint(fingerprints_stack[-1], child_index) if fingerprints_stack
                    else fingerprints_)

            # len(path_interval_sums) > MAX_CONTINUE_FN + 1:
            if fn_factors[-1] < CRITICAL_FN_FACTOR:
                continue
            for child_index, (child_site, interval, _, _) in enumerate(current_site.children):
                debug = False
                result.add(child_site)
                frontier.append((child_index, interval, child_site))
                whole_path_length = path_interval_sums[-1] + interval

                prob_candidate = content * (1 - prob_skip(whole_path_length)) * fn_factors[-1] * similar_factor(whole_path_length, intervals_)
                fingerprints_candidate = update_fingerprint(fingerprints_stack[-1], child_index)
                print('{} --{},{}-> {}'.format(start_site.id, content.shape[0], len(child_indexs), child_site.id))
                modify(prob_matrixs[child_site][i], tracker_matrixs[child_site][i], fingerprint_matrixs[child_site][i],
                        prob_candidate, fingerprints_candidate, start_site, child_indexs[1:] + [child_index])
        else:
            for ele in site_path, child_indexs, fn_factors, path_interval_sums, fingerprints_stack:
                ele.pop()
    return result

def trace_back(prob_matrixs, tracker_matrixs, mol_index, tracker):
    """
    Trace back a route start from a tracker.
    Arguments:
        prob_matrixs (numpy.array, size: m x n): The coorespond prob table.
        tracker_matrixs (numpy.array, size: m x n): The coorespond tracker table.
        mol_index (int): the index of row which the tracker is in.
        tracker (tuple): initial tracker.
    Return:
        (tuple)
    """
    actions, site_path, node_path, probs, lengths, child_indexs = [], [], [], [], [], []
    print('tracker:', tracker)
    print('mol index:', mol_index)
    while tracker:
        pre_site, (delta_mol_index, pre_tracker_rank), sub_child_indexs = tracker
        print('tracker:', tracker)
        print('pre_site:', pre_site)
        print(delta_mol_index, pre_tracker_rank)
        print('sub_child_indexs:', sub_child_indexs)
        for i in range(len(sub_child_indexs) - 1, -1, -1):
            actions.append('{}-{}'.format('M' if i == 0 else 'MM', sub_child_indexs[i]))
        actions.extend(['S'] * (-delta_mol_index - 1))
        mol_index += delta_mol_index
        probs.append(prob_matrixs[pre_site][mol_index][pre_tracker_rank])

        length = 0
        sub_site_path, sub_node_path = [pre_site], []
        current_site = pre_site
        for child_index in sub_child_indexs:
            child_site, interval, sub_sub_node_path, contamination = current_site.children[child_index]
            sub_sub_node_path = list(zip(sub_sub_node_path, contamination))
            if sub_node_path:
                sub_node_path.extend(sub_sub_node_path[1:])
            else:
                sub_node_path.extend(sub_sub_node_path)
            print(interval, 'nodes:', [(ele[0].uid, ele[1]) for ele in sub_sub_node_path])
            sub_site_path.append(child_site)
            length += interval
            current_site = child_site
        if site_path:
            lengths.append(length)
            sub_site_path.pop()
        if node_path:
            sub_node_path.pop()
        sub_site_path.reverse()
        sub_node_path.reverse()
        site_path.extend(sub_site_path)
        node_path.extend(sub_node_path)
        child_indexs.extend(sub_child_indexs[::-1])
        tracker = tracker_matrixs[pre_site][mol_index][pre_tracker_rank]
        prob = prob_matrixs[pre_site][mol_index][pre_tracker_rank]
        print('mol_index:', mol_index)
        print('prob of last tracker:', prob)
        print('tracker:', tracker)
    print('mol index:', mol_index)
    print()
    assert mol_index == 0
    for list_ in actions, site_path, probs, node_path, lengths, child_indexs:
        list_.reverse()
    # print('hello')
    # print('node_path:', node_path)
    print('node_path:', node_path)
    node_path, contamination = zip(*node_path)
    return site_path, node_path, contamination, probs, lengths, actions, child_indexs

def find_path(start_sites, end_sites, sites, intervals, n_rank):
    """
    Find `n_rank` paths in site graph that the interval of sites on path best match `intervals`.
    Argument:
        start_sites (list[site_graph.Site]): Possible start sites of the wanted path.
        end_sites (list[site_graph.Site]): Possible end sites of the wanted path.
        sites (dict[str:site_graph.Site]): site graph.
        intervals (list[float]): intervals we want the path to match.
        n_rank (int): number of path to output.
    Note:
        When `start_sites` or `end_sites` is empty, it means all sites included.
    Return:
        (list[tuple])
    """
    if not start_sites:
        start_sites = sites.values()
    if not end_sites:
        end_sites = sites.values()
    result = []
    len_intervals = len(intervals)
    # Create one pro_matrix and one tracker_matrix for each site.
    prob_matrixs, tracker_matrixs, fingerprint_matrixs = {}, {}, {}
    matrix_size = (len_intervals + 1, n_rank)
    for site in sites.values():
        prob_matrixs[site] = np.zeros(matrix_size)
        tracker_matrixs[site] = np.empty(matrix_size, dtype='object')
        fingerprint_matrixs[site] = np.zeros(matrix_size, dtype='uint64')
    start_prob = 1 / len(start_sites)
    for start_site in start_sites:
        prob_matrixs[start_site][0][0] = start_prob
    fingerprint_matrixs[start_site][0][0] = hash(start_site.id)
    reached_sites = set(start_sites)
    accu_product_log10 = 0
    for i in range(1, matrix_size[0]):
        print('i:', i)
        new_reached_sites = set()
        print("len_reached_sites:", len(reached_sites))
        for site in reached_sites:
            new_reached_sites.update(propagate(i, site, intervals, prob_matrixs, tracker_matrixs, fingerprint_matrixs, 1 / len(sites) * 1e-5))
            # new_reached_sites.update(propagate(i, site, intervals, prob_matrixs, tracker_matrixs, fingerprint_matrixs, 0))
        # print('new reached sites')
        # print([ele.id for ele in new_reached_sites])
        reached_sites.update(new_reached_sites)
        print('number of reached sites:', len(reached_sites))
        print([ele.id for ele in reached_sites])
        # Normalize the ith row of each prob_matrix.
        total_sum = sum((ele[i][0] for ele in prob_matrixs.values()))
        print('total_sum:', total_sum)
        accu_product_log10 += math.log10(total_sum)
        total_sum_invert = 1 / total_sum
        for ele in prob_matrixs.values():
            ele[i] *= total_sum_invert

    # f_prob = open('prob_matrixs_0', 'wb')
    # f_tracker = open('tracker_matrixs_0', 'wb')
    
    # prob_matrixs_pickle, tracker_matrixs_pickle = {}, {}
    # for key, value in prob_matrixs.items():
    #     prob_matrixs_pickle[key.id] = value
    # for key, value in tracker_matrixs.items():
    #     new_matrix = np.empty_like(value)
    #     for i in range(value.shape[0] * value.shape[1]):
    #         old_content = value.flat[i]
    #         new_content = (old_content[0].id, old_content[1], old_content[2]) if old_content else None
    #         new_matrix.flat[i] = new_content
    #     tracker_matrixs_pickle[key.id] = new_matrix
    # pickle.dump(prob_matrixs_pickle, f_prob, -1)
    # pickle.dump(tracker_matrixs_pickle, f_tracker, -1)
 
    # f_prob.close()
    # f_tracker.close()

    site_path_probs, final_trackers, final_fingerprints = np.zeros(n_rank), np.empty(n_rank, dtype='object'), np.empty(n_rank, dtype='uint64')
    for end_site in end_sites:
        modify(site_path_probs, final_trackers, final_fingerprints,
                # get_latest_row(prob_matrixs[end_site], matrix_size[0]), get_latest_row(fingerprint_matrixs[end_site], matrix_size[0]), end_site, [])
                prob_matrixs[end_site][matrix_size[0]-1:matrix_size[0], :], fingerprint_matrixs[end_site][matrix_size[0]-1:matrix_size[0], :], end_site, [])

    print('site_path_probs:', site_path_probs)
    print('final_tracker:', final_trackers)
    print('sum:', site_path_probs.sum())
    accu_product_log10 += math.log10(site_path_probs.sum())
    site_path_probs = site_path_probs / site_path_probs.sum()
    print('site_path_probs:', site_path_probs)
    # print('accu_product, len_interavl, log/len:', accu_product, len(intervals), math.log10(accu_product) / len(intervals))
    print('accu_product_log10, len_interavl, log/len:', accu_product_log10)
    print(len(intervals))
    print(accu_product_log10 / len(intervals))

    total_prob = .999
    index = 0
    accu_prob = site_path_probs[0]
    while accu_prob < total_prob:
        index += 1
        accu_prob += site_path_probs[index]

    # Trace back.
    for path_prob, tracker in zip(site_path_probs[:index + 1], final_trackers[:index + 1]):
        result.append((path_prob,) + trace_back(prob_matrixs, tracker_matrixs, matrix_size[0], tracker))
    return result

def get_measure_len(action_strs, intervals):
    result = []
    action_strs = [ele for ele in action_strs if 'MM' not in ele]
    reserved_len = 0
    for action_str, interval in zip(action_strs, intervals):
        if 'M' in action_str:
            result.append(interval + reserved_len)
            reserved_len = 0
        else:
            reserved_len += interval
    return result

def get_seq(node_path, is_valids, start_trim_len, end_keep_len):
    is_valids = [True] * len(is_valids)
    overlap = start_trim_len if start_trim_len else 0
    current_node = node_path[0]
    is_last_node_valid = True
    result = []
    for i in range(len(node_path)):
        assert current_node == node_path[i]
        # result.append(current_node.seq if is_valids[i] else 'N' * len(node_path[i].seq))
        if not is_last_node_valid:
            # Shrink `overlap` Ns from last N seq.
            result.append(result.pop()[overlap:])
            overlap = 0
        result.append(current_node.seq[overlap:] if is_valids[i] else (len(current_node.seq) - overlap) * 'N')
        if i < len(node_path) - 1:
            j = 0
            while node_path[i + 1] != current_node.children[j][0]:
                j += 1
            current_node, overlap = current_node.children[j]
            is_last_node_valid = is_valids[i]
    if end_keep_len:
        result.append(result.pop()[:end_keep_len - overlap])
    return ''.join(result)

def process_find_path_result(result, sites, intervals):
    node_paths = []
    accu_prob = 0
    i = 0
    contaminations = []
    for path_prob, site_path, node_path, contamination, probs, lengths, actions, child_indexs in result:
        if accu_prob > .95:
            break
        print('i:', i)
        i += 1
        measured_lens = get_measure_len(actions, intervals)
        print('actual len:', lengths)
        print('measured_lens:', measured_lens)
        actual_num_FN = sum([(1 if 'MM' in ele else 0) for ele in actions])
        print('actual FN num:', actual_num_FN)
        actual_num_FP = sum([(1 if 'S' in ele else 0) for ele in actions])
        print('expected:', round(sum([(1 if 'M' in ele else 0) for ele in actions]) * FN_RATE, 2))
        print('actual FP num:', actual_num_FP)
        print('expected:', round(sum(intervals) / AVERAGE_INTERVAL, 2))
        print('site_path:', len(site_path), [ele.id for ele in site_path])
        print('child_indexs:', len(child_indexs), child_indexs)
        errors = [(measured_len - actual_len) for (measured_len, actual_len) in zip(measured_lens, lengths)]
        print('errors:', errors)
        average = statistics.mean(errors)
        if len(errors) < 2:
            std = abs(errors[0] - average)
        else:
            std = statistics.stdev(errors, average)
        print('actual average & std:', average, std)
        print('expected:', MU, SIGMA)
        accu_prob += path_prob
        contaminations.append(contamination)
        node_paths.append((node_path, contamination))
        print('node_path:', [(ele.uid, contaminated) for (ele, contaminated) in zip(node_path, contamination)])
    print('accu_prob:', accu_prob)
    is_valids = merge_node_path(node_paths)
    for node, is_valid in zip(node_paths[0][0], is_valids):
        print(node, is_valid)
    return node_paths[0][0], is_valids


def positions_to_intervals(positions):
    return [(positions[i] - positions[i - 1]) for i in range(1, len(positions))]

def merge_node_path(paths):
    paths_len = len(paths)
    template_path, template_contamination = paths[0]
    matcher = difflib.SequenceMatcher()
    matcher.set_seq1(template_path)
    support = [0 if contaminated else 1 for contaminated in template_contamination]
    for path, contamination in paths[1:]:
        matcher.set_seq2(path)
        print(matcher.get_matching_blocks())
        for index, index_q, length in matcher.get_matching_blocks():
            while length > 0:
                if not(template_contamination[index] or contamination[index_q]):
                    support[index] += 1
                index += 1
                index_q += 1
                length -= 1
    return [ele == paths_len for ele in support]

def print_help():
    pass

def main():
    n_rank = 10
    options, args = getopt.getopt(sys.argv[1:], 'c:s:e:a:b:i:g:o:r:n:')
    start_site_ids, end_site_ids = [], []
    position_in_start_node, position_in_end_node = None, None
    gap_name = None
    for option, value in options:
        if option == '-c':
            MU, SIGMA = map(float, value.split(','))
        elif option == '-s':
            start_site_ids = value.split(',')
        elif option == '-e':
            end_site_ids = value.split(',')
        elif option == '-n':
            gap_name = value
        elif option == '-a':
            position_in_start_node = int(float(value))
        elif option == '-b':
            position_in_end_node = int(float(value))
        elif option == '-i':
            measures = list(map(float, value.split(',')))
        elif option == '-g':
            graph_file = value
            with open(graph_file, 'rb') as fin:
                nodes = pickle.load(fin)
                sites = pickle.load(fin)
        elif option == '-o':
            output_file_name = value
        elif option == '-r':
            n_rank = int(value)

    start_sites = [sites[ele] for ele in start_site_ids]
    end_sites = [sites[ele] for ele in end_site_ids]
    result = find_path(start_sites, end_sites, sites, measures, n_rank)
    node_path, is_valids = process_find_path_result(result, sites, measures)
    with open(output_file_name, 'a') as fout:
        fout.write('>{}\n'.format(gap_name))
        fout.write(get_seq(node_path, is_valids, position_in_start_node, position_in_end_node))
        fout.write('\n')


if __name__ == '__main__':
    tick = time.time()
    main()
    # _test_get_seq()
    # fastg_file_name, site_graph_file_name = sys.argv[1:3]
    # _test_find_path3(fastg_file_name, site_graph_file_name, overlap=77)
    # _test_merge_node_path()
    # _test_similar_factor()
    # _test_prob_skip()
    tock = time.time()
    print('Time used:', round(tock - tick, 2), 'seconds')


