"""
This script is a rewrite of file find_path_dp.py, aims to speedup, better logging,
and clearer output.
"""

import logging
import argparse
import sys
import pickle
import numpy as np
from numba import njit

# Global variable.
logger = None
## Default parameter.
MU = 0
SIGMA = 300
RANK = 10

HASH_A = 4343534
HASH_B = 232423
def update_fingerprint(original_fingerprint, child_index):
    return (original_fingerprint + HASH_A) * (child_index + HASH_B)

class BackTracer:

    def __init__(self, site_id='', x=-1, y=-1, child_index=-1):
        self.site_id = site_id
        self.x = x
        self.y = y
        self.child_index = child_index

    @property
    def is_start(self):
        return True if self.site_id == '' else False

class PathFinder:

    def __init__(self, nrank):
        self._n_rank = nrank
        self._p_tensor = None
        self._t_tensor = None
        self._f_tensor = None
        self._start_sites = set()
        self._end_sites = set()
        self._site_id_to_index = {}

    def load_graph_and_interval(self, sites, intervals):
        logger.info("Loading site graph and intervals ...")
        shape = (len(sites), len(intervals) + 1, self._n_rank)
        self._p_tensor = np.zeros(shape, dtype='float32')
        self._t_tensor = np.empty(shape, dtype='object')
        self._f_tensor = np.zeros(shape, dtype='uint64')
        logger.info("Loaded %d sites, %d intervals.", len(sites), len(intervals))
        site_ids = sorted(list(sites.keys()))
        self._site_id_to_index = {k: v for v, k in enumerate(site_ids)}

    def site_p_tensor(self, site_id):
        return self._p_tensor(self._site_id_to_index[site_id])

    def site_t_tensor(self, site_id):
        return self._t_tensor(self._site_id_to_index[site_id])

    def site_f_tensor(self, site_id):
        return self._f_tensor(self._site_id_to_index[site_id])

    def loading_start_end_sites(self, start_sites=[], end_sites=[]):
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
    
    def run(self):
        pass


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

    logger.info('Reading graph files ...')
    nodes = pickle.load(args.graph_file)
    sites = pickle.load(args.graph_file)
    logger.info('Loaded an assembly graph, containing %d nodes.', len(nodes))
    logger.info('Loaded a site graph, containing %d sites.', len(sites))
    # args.graph_file.close()  # Needed?

    finder = PathFinder(args.rank)
    finder.load_graph_and_interval(sites, args.intervals)
    finder.loading_start_end_sites(args.start_sites, args.end_sites)
    

if __name__ == "__main__":
    main()