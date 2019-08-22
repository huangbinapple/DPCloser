"""
This script is a rewrite of file find_path_dp.py, aims to speedup, better logging,
and clearer output.
"""

import logging
import argparse
import sys
import math
import pickle
import numpy as np
from numba import njit, vectorize

# Global variable.
logger = None
## Default parameter.
MU = 0
SIGMA = 300
RANK = 10
MERGE_LEN = 500
FN_RATE = .228
MAX_CONTINUE_FN = 1.9
MAX_CONTINUE_FP = 2
AVERAGE_INTERVAL = 1526000

FIVE_SIGMA = 5 * SIGMA
DECAY_INVERSE = math.sqrt(2) / SIGMA
CRITICAL_FN_FACTOR = FN_RATE ** (MAX_CONTINUE_FN)

HASH_A = 4343534
HASH_B = 232423
def poisson(lambda_, k):
    return math.exp(-lambda_) * lambda_ ** k / math.factorial(k)

def find_none(ndarray):
    size = ndarray.shape[0]
    i = 0
    while i < size:
        if None is ndarray[i]:
            return i
        i += 1
    return i

class BackTracer:

    def __init__(self, site_index=-1, delta_x=0, delta_y=-1, child_indexs=[]):
        self.site_index = site_index
        self.delta_x = delta_x
        self.delta_y = delta_y
        self.child_indexs = child_indexs

    @property
    def is_start(self):
        return True if self.site_index == -1 else False

    def __str__(self):
        return str((self.site_index, self.delta_x, self.delta_y, self.child_indexs))

    def __repr__(self):
        return str(self)

class PathFinder:

    def __init__(self, nrank):
        self._n_rank = nrank
        self._p_tensor = None
        self._t_tensor = None
        self._f_tensor = None
        self._intervals = None
        self._sites = None
        self._start_sites = set()
        self._end_sites = set()
        self._site_id_to_index = {}
        self._site_ids = []
        self._log_sum = 0
        self._child_index = {}
        self._propagation_index = {}
        self._interval_index = []

    def load_graph(self, sites):
        logger.info("Loading site graph...")
        self._sites = sites
        logger.info("Loaded %d sites.", len(sites))
        self._site_ids = sorted(list(sites.keys()))
        self._site_id_to_index = {k: v for v, k in enumerate(self._site_ids)}
        self.index_graph()

    def load_intervals(self, intervals):
        assert self._sites
        self._intervals = intervals
        logger.info("Loading site intervals...")
        shape = (len(self._sites), len(intervals) + 1, self._n_rank)
        self._p_tensor = np.zeros(shape)
        self._t_tensor = np.empty(shape, dtype='object')
        self._f_tensor = np.zeros(shape, dtype='uint64')
        logger.info("Loaded %d intervals.", len(intervals))
        self.index_intervals()

    def site_ids(self, index):
        return self._site_ids[index]

    def site_p_tensor(self, site_id):
        return self._p_tensor(self._site_id_to_index[site_id])

    def site_t_tensor(self, site_id):
        return self._t_tensor(self._site_id_to_index[site_id])

    def site_f_tensor(self, site_id):
        return self._f_tensor(self._site_id_to_index[site_id])
    
    @staticmethod
    def similar_factor(length_reference, length_bionano):
        """
        Argument:
            length_reference (numpy.array, size: m x 1, dtype: int):
                path length in the site graph.
            length_bionano (numpy.array, size: 1 x n):
                intervals on bionano molecules. 
        Return:
            (numpy.array, size: `length_reference` operate `length_bionano`):
                score factor represent the similarity between measured and reference length.
        """
        error = np.abs(length_bionano - (MU + length_reference))
        # error = np.minimum(error, FIVE_SIGMA)
        result = np.exp(- error * DECAY_INVERSE)
        return result

    @staticmethod
    @vectorize
    def prob_skip(interval):
        if interval < MERGE_LEN:
            return (1 - FN_RATE) * (interval / MERGE_LEN - 1) ** 2 + FN_RATE
        else:
            return FN_RATE

    def _index_children(self, site):
        site_index = self._site_id_to_index[site.id]
        if site.children:
            children_zip = zip(*site.children)
            # children = np.array(next(children_zip))
            children = np.array(list(map(self._site_id_to_index.get, map(lambda x: x.id, next(children_zip)))))
            intervals = np.array(next(children_zip))
            probs_skip = self.prob_skip(intervals)
            self._child_index[site_index] = children, intervals, probs_skip
        else:
            self._child_index[site_index] = (
                np.empty(0, dtype='int64'),
                np.empty(0, dtype='int64'),
                np.empty(0, dtype='float64')
            )

    def _index_propagation_route(self, site):
        logger.debug('Building propagation route for site %s.', site.id)
        sub_result = []
        n_iter = 0
        index = self._site_id_to_index[site.id] 
        # Init loop variables.
        indexs, intervals, _ = self._child_index[index]
        # sites, intervals, _ = self._child_index[site]
        init_size = len(indexs)
        children_indexs = np.empty(init_size, dtype=object)
        for i in range(len(children_indexs)):
            children_indexs[i] = [i]
        proceed_lower_bounds = np.full(init_size, CRITICAL_FN_FACTOR)
        logger.debug('Init. %d edges.', init_size)
        while True:
            to_proceeds = np.array([proceed_lower_bounds[i] < self._child_index[indexs[i]][2]
                for i in range(len(indexs))])  # True means be able to propagate.
            num_to_proceeds = np.array([np.count_nonzero(ele) for ele in to_proceeds])
            index_propagate = (num_to_proceeds > 0)
            # Dump.
            fn_scores = CRITICAL_FN_FACTOR / proceed_lower_bounds
            sub_result.append((indexs, children_indexs, intervals, fn_scores))
            new_size = num_to_proceeds.sum()
            logger.debug('#Length %d edges dumped: %d.', n_iter + 1, np.count_nonzero(indexs))
            if new_size == 0:
                break
            ## Propagate.
            logger.debug('Propagating from %d edges(%d edges when finished.)...',
                np.count_nonzero(index_propagate), num_to_proceeds.sum())
            # Construct loop variable for next loop.
            new_indexs = np.empty(new_size, dtype='int64')
            new_intervals = np.empty(new_size, dtype='int64') 
            new_children_indexs = np.empty(new_size, dtype=object)
            new_proceed_lower_bounds = np.empty(new_size, dtype='float64')
            start_index = 0
            for i in index_propagate.nonzero()[0]:
                to_proceed = to_proceeds[i]
                end_index = start_index + num_to_proceeds[i]
                next_indexs, next_intervals, next_probs_skip = self._child_index[indexs[i]]
                new_indexs[start_index:end_index] = next_indexs[to_proceed]
                new_intervals[start_index:end_index] = intervals[i] + next_intervals[to_proceed]
                for delta, next_child_index in enumerate(to_proceed.nonzero()[0]):
                    new_children_indexs[start_index + delta] = children_indexs[i] + [next_child_index]
                new_proceed_lower_bounds[start_index:end_index] = proceed_lower_bounds[i] / next_probs_skip[to_proceed]
                start_index = end_index
            # Update loop variable.
            indexs, intervals, children_indexs, proceed_lower_bounds =\
                new_indexs, new_intervals, new_children_indexs, new_proceed_lower_bounds
            n_iter += 1
        self._propagation_index[index] = tuple((np.concatenate(ele) for ele in zip(*sub_result)))

    @staticmethod
    # @njit
    def _normalize(index, tensor):
        logger.debug('Normalize on: %s', tensor[:, index, 0])
        logger.debug('Normalize on(sorted): %s', np.sort(tensor[:, index, 0])[::-1])
        logger.debug("All zeros: %s", not any(tensor[:, index, 0]))
        sum_ = tensor[:, index, 0].sum()
        sum_inverse = 1 / sum_
        tensor[:, index, :] *= sum_inverse
        logger.debug('Row %d sum: %d', index, sum_)
        logger.debug('Row %d sum_inverse: %d', index, sum_inverse)
        return sum_
    
    def normalize(self, index):
        total_sum = self._normalize(index, self._p_tensor)
        self._log_sum += np.log2(total_sum)
        return total_sum

    @staticmethod
    def propagate_fingerprints(init_fingerprints, children_indexs):
        logger.debug('Updating fingerprnits: %s ...', init_fingerprints)
        result_fingerprints = np.empty(
            shape=(children_indexs.shape[0], init_fingerprints.shape[0]),
            dtype='uint64'
        )
        last_children_index_len = 0
        memory = {}  # Values store the index of updated fingerprints in result_fingerprints.
        for i, children_indexs_ in enumerate(children_indexs):
            current_children_index_len = len(children_indexs_)
            assert current_children_index_len >= last_children_index_len
            last_children_index_len = current_children_index_len
            if current_children_index_len > 1:
                init_fingerprints_ = result_fingerprints[memory[tuple(children_indexs_[:-1])], :]
            else:
                init_fingerprints_ = init_fingerprints
            memory[tuple(children_indexs_)] = i
            result_fingerprints[i] = PathFinder.update_fingerprint(
                init_fingerprints_, children_indexs_[-1])
        logger.debug('%d fingerprint update completed!', len(children_indexs))
        return result_fingerprints

    @staticmethod
    @njit
    def update_fingerprint(init_fingerprints, child_index):
        return (init_fingerprints + HASH_A) * (child_index + HASH_B)
    
    def find_path(self):
        logger.info('Finding optimal paths ...')
        log_sum = self._find_path(self._p_tensor, self._t_tensor, self._f_tensor,
            self._site_ids, self._interval_index, self._propagation_index, mini_prob=1e-5)
        self._log_sum += log_sum
        logger.info('Find paths complete, log sum is %d', self._log_sum)

    @staticmethod
    # @njit
    def _find_path(P, T, F, site_ids, interval_index, propagate_index, mini_prob=0):
        index_iter = 0
        total_iter = len(interval_index)
        log_sum = 0
        logger.info('Probability distribution of #%d site in site graph nodes %s', 0,
            sorted(list(zip([site_ids[id_] for id_ in P[:, 0, 0].nonzero()[0]],
            P[:, index_iter, 0][P[:, 0, 0].nonzero()[0]])), key=lambda x: x[1], reverse=True))
        for index_iter in range(total_iter):
            logger.info('Finding progress %d/%d.', index_iter + 1, total_iter)
            chosen_index = (P[:,index_iter,0] > mini_prob).nonzero()[0]
            logger.info('%d site to propagate.', len(chosen_index))
            for site_index in chosen_index:
                PathFinder.propagate(P, T, F, site_ids, 
                    site_index, index_iter, interval_index[index_iter], propagate_index[site_index])
            normalize_sum = PathFinder._normalize(index_iter + 1, P)
            logger.info('Probability distribution of #%d site in site graph nodes %s', index_iter + 1,
                sorted(list(zip([site_ids[id_] for id_ in P[:, index_iter + 1, 0].nonzero()[0]],
                P[:, index_iter + 1, 0][P[:, index_iter + 1, 0].nonzero()[0]])), key=lambda x: x[1], reverse=True))
            log_sum += np.log2(normalize_sum)
            # break
        return log_sum
            
    @staticmethod
    # @njit
    def propagate(P, T, F, site_ids, site_index, index_iter, interval_index_content, propogate_index_content):
        # Input.
        site_id = site_ids[site_index]
        logger.info('(Iter: %d) Propagating from site %s...', index_iter, site_id)
        num_propagate = find_none(T[site_index][index_iter])
        logger.debug("num_propagate: %d", num_propagate)
        init_probs = P[site_index][index_iter]
        init_fingerprints = F[site_index][index_iter]
        target_indexs, children_indexs, reference_lengths, fn_scores = propogate_index_content
        bionano_lengths, fp_scores = interval_index_content

        propagate_factor = PathFinder.similar_factor(reference_lengths.reshape(-1, 1), bionano_lengths) * \
            fn_scores.reshape(-1, 1) * fp_scores
        logger.debug('Defactor propagate factor:')
        logger.debug('reference lengths: %s', reference_lengths.reshape(-1, 1))
        logger.debug('Bionano lengths: %s', bionano_lengths)
        logger.debug('Similate factor: %s', PathFinder.similar_factor(reference_lengths.reshape(-1, 1), bionano_lengths))
        logger.debug('fn_scores: %s', fn_scores.reshape(-1, 1))
        logger.debug('fp_scores: %s', fp_scores)
        logger.debug('Propagate factors: %s', propagate_factor)
        propagate_results = np.expand_dims(propagate_factor, axis=-1) * init_probs
        updated_fingerprints = PathFinder.propagate_fingerprints(init_fingerprints, children_indexs)

        # Action(Alter table.)
        for target_index, children_indexs_, propagate_result, fingerprints in zip(
                target_indexs, children_indexs, propagate_results, updated_fingerprints):
            tracker = T[target_index][index_iter + 1]
            num_already_here = find_none(tracker)
            logger.debug("Alter table operation, target site id %s:", site_ids[target_index])
            logger.debug('num_already_here: %d', num_already_here)
            logger.debug("Original P: %s", P[target_index][index_iter + 1])
            logger.debug("Original F: %s", F[target_index][index_iter + 1])
            logger.debug("Original T: %s", T[target_index][index_iter + 1])
            logger.debug("Children indexs: %s", children_indexs_)
            logger.debug("Propagate result: %s", propagate_result)
            logger.debug("Fingerprint: %s", fingerprints)
            tracker_info = PathFinder.merge(P[target_index][index_iter + 1], propagate_result,
                F[target_index][index_iter + 1], fingerprints, num_already_here, num_propagate)
            # Update tracker.
            tracker = tracker[:tracker_info.shape[0]]
            tracker_keep_index = (tracker_info < num_already_here)
            tracker[tracker_keep_index] = tracker[tracker_info[tracker_keep_index]]
            tracker[~ tracker_keep_index] = tuple(map(
                lambda x:
                BackTracer(
                    site_index,
                    - (x - num_already_here) // num_propagate,
                    (x - num_already_here) % num_propagate,
                    children_indexs_
                ),
                tracker_info[~ tracker_keep_index]
            ))
            logger.debug("Finish a alter table operation, tracker info: %s", tracker_info)
            logger.debug("Altered P: %s", P[target_index][index_iter + 1])
            logger.debug("Altered F: %s", F[target_index][index_iter + 1])
            logger.debug("Altered T: %s", T[target_index][index_iter + 1])
             
        logger.debug('Propagation finished!')

    @staticmethod
    @njit
    def merge(values_a, values_b, fingerprints_a, fingerprints_b, n_valid_a, n_valid_b):
        size = values_a.shape[0]
        # logger.debug('n_valid_a: %d', n_valid_a)
        # logger.debug('n_valid_b: %d', n_valid_b)
        assert len(values_b.shape) == 2
        b_row = values_b.shape[0]
        values_ab = np.concatenate((values_a[:n_valid_a], values_b[:, :n_valid_b].flatten()))
        # logger.debug("value_ab: %s", values_ab)
        fingerprints_b_repeat = np.empty((b_row, n_valid_b), dtype=np.uint64)
        for i in range(b_row):
            fingerprints_b_repeat[i] = fingerprints_b[:n_valid_b]
        fingerprints_ab = np.concatenate((fingerprints_a[:n_valid_a], fingerprints_b_repeat.flatten()))
        # logger.debug("fingerprints_ab: %s", fingerprints_ab)
        index_sort = (-values_ab).argsort(kind='mergesort')  # Descend order.
        # logger.debug("index_sort: %s", index_sort)
        fingerprints_sorted_by_value = fingerprints_ab[index_sort]
        # logger.debug("fingerprints_sorted_by_value: %s", fingerprints_sorted_by_value)

        fingerprint_sort_index = fingerprints_sorted_by_value.argsort(kind='mergesort')  # Merge sort matters
        # logger.debug("fingerprint_sort_index: %s", fingerprint_sort_index)
        fingerprint_sorted = fingerprints_sorted_by_value[fingerprint_sort_index]
        # logger.debug("fingerprint_sorted: %s", fingerprint_sorted)
        flag = np.empty(len(fingerprint_sort_index), dtype=np.bool_)
        flag[0] = True
        flag[1:] = (fingerprint_sorted[1:] != fingerprint_sorted[:-1])
        index_unique = fingerprint_sort_index[flag]
        index_unique.sort()
        # logger.debug("index_unique: %s", index_unique)

        index_result = index_sort[index_unique[:size]]
        # logger.debug("index_result: %s", index_result)
        result_size = index_result.shape[0]
        # logger.debug('result size: %d', result_size)
        # logger.debug("Values altered in: %s", values_ab[index_result])
        values_a[:result_size] = values_ab[index_result]
        # logger.debug("Value alterd in data type: %s", values_ab.dtype)
        # logger.debug("Value_a: %s", values_a)
        # logger.debug("All zeros: %s", not any(values_a))
        # logger.debug("Value a data type: %s", values_a.dtype)
        fingerprints_a[:result_size] = fingerprints_ab[index_result]
        # logger.debug("Fingerprint altered in: %s", fingerprints_ab[index_result])
        return index_result

    def load_start_end_sites(self, start_sites=[], end_sites=[]):
        assert not self._start_sites and not self._end_sites
        if not start_sites:
            start_sites = set(self._site_id_to_index.values())
        if not end_sites:
            end_sites = set(self._site_id_to_index.values())
            
        logger.info('Loading start sites and end sites: %s, %s.', str(start_sites), str(end_sites))
        self._start_sites.update(start_sites)
        self._end_sites.update(end_sites)
        init_prob = 1 / len(start_sites)
        for site_index in map(self._site_id_to_index.get, start_sites):
            self._p_tensor[site_index][0][0] = init_prob
            self._t_tensor[site_index][0][0] = BackTracer()
            self._f_tensor[site_index][0][0] = site_index
        logger.info('Loaded %d start sites and %d end sites.', len(start_sites), len(end_sites))
        self.normalize(0)
    
    def index_graph(self):
        logger.info('Indexing graph...')
        logger.info('Build child index for all sites...')
        for site in self._sites.values():
            self._index_children(site)
        logger.info('All site child index built.')
            
        logger.info('Build propagate for all sites...')
        for site in self._sites.values():
            self._index_propagation_route(site)
        logger.info('All site propagation index built.')

        count = sum((len(ele[0]) for ele in self._propagation_index.values()))
        logger.info('Graph index built, %d site, %d edges.', len(self._propagation_index), count)

    def index_intervals(self):
        logger.info('Indexing intervals...')
        intervals_ = self._intervals[::-1]
        result = []
        for i in range(len(intervals_)):
            interval = 0
            intervals = []
            fp_scores = []
            for num_insert in range(MAX_CONTINUE_FP):
                if i + num_insert >= len(intervals_):
                    break
                interval += intervals_[i + num_insert]
                intervals.append(interval)
                fp_score = poisson(interval / AVERAGE_INTERVAL, num_insert)
                fp_scores.append(fp_score)
                logger.debug('Prob of FN=%d in %d is %f', num_insert, interval, fp_score)
            result.append((np.array(intervals), np.array(fp_scores)))
            self._interval_index = result[::-1]
        logger.info('Index intervals done.')
        logger.debug('Interval index: %s', self._interval_index)

def main():
    global MU
    global SIGMA
    global RANK
    global logger

    parser = argparse.ArgumentParser()

    parser.add_argument('graph_file',
        type=argparse.FileType('rb'),
        metavar='INPUT_FILE', help='A pickle file.')

    parser.add_argument('--start_sites',
        required=True, type=lambda x: [ele for ele in x.split(',')],
        metavar='START_SITE_ID', help='All possible start site ids, split by comma.')

    parser.add_argument('--end_sites',
        required=True, type=lambda x: [ele for ele in x.split(',')],
        metavar='END_SITE_ID', help='All possible end site ids, split by comma.')

    parser.add_argument('--intervals',
        required=True, type=lambda x: [int(ele) for ele in x.split(',')],
        metavar='INTERVAL', help='All intervals between two neighboring sites, split by comma.')

    parser.add_argument('--mu',
        default=MU, type=int,
        metavar='MU', help='Expected shift in measurement relative to ref.')

    parser.add_argument('--sigma',
        default=SIGMA, type=int,
        metavar='SIGMA', help='Variance of shift in measurement relative to ref.')

    parser.add_argument('--rank',
        default=RANK, type=int,
        metavar='RANK', help='Variance of shift in measurement relative to ref.')

    parser.add_argument('gap_name',
        metavar='GAP_NAME', help='The gaps name, give a name whatever you want.')

    parser.add_argument("-l", "--log", dest="log_level",
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='INFO',
        metavar='LOG_LEVEL', help="Set the logging level.")

    args=parser.parse_args()

    MU = args.mu
    SIGMA = args.sigma
    RANK = args.rank

    # Set logger.
    logger = logging.getLogger('path_finder')
    logger.setLevel(getattr(logging, args.log_level))
    file_stream = logging.FileHandler(args.gap_name + '.log', mode='w')
    file_stream.setLevel(getattr(logging, args.log_level))
    file_stream.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    ))
    logger.addHandler(file_stream)
    
    logger.info('Input command: %s', ' '.join(sys.argv))
    logger.info('Gap name: %s', args.gap_name)
    logger.info('SIGMA: %d', SIGMA)
    logger.info('RANK: %d', RANK)
    logger.info('Start site ids (%d ids): %s', len(args.start_sites), args.start_sites)
    logger.info('End site ids (%d ids): %s', len(args.end_sites), args.end_sites)
    logger.info('Measurements (%d intervals): %s', len(args.intervals), args.intervals)

    logger.info('Reading graph files...')
    nodes = pickle.load(args.graph_file)
    sites = pickle.load(args.graph_file)
    logger.info('Loaded an assembly graph, containing %d nodes.', len(nodes))
    logger.info('Loaded a site graph, containing %d sites.', len(sites))
    # args.graph_file.close()  # Needed?

    finder = PathFinder(args.rank)
    finder.load_graph(sites)
    finder.load_intervals(args.intervals)
    finder.load_start_end_sites(args.start_sites, args.end_sites)
    finder.find_path()


if __name__ == "__main__":
    main()