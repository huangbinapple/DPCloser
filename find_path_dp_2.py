"""
This script is a rewrite of file find_path_dp.py, aims to speedup, better logging,
and clearer output.
"""

import logging
import argparse
import sys
import pickle

# Default parameter.
MU = 0
SIGMA = 300

def main():
    global MU
    global SIGMA

    parser = argparse.ArgumentParser()

    parser.add_argument('graph_file',
        metavar='INPUT_FILE', help='A pickle file.')

    parser.add_argument('--start_sites',
        required=True, nargs='+',
        metavar='START_SITE_ID', help='All possible start site ids.')

    parser.add_argument('--end_sites',
        required=True, nargs='+',
        metavar='END_SITE_ID', help='All possible end site ids.')

    parser.add_argument('--intervals',
        required=True, type=float, nargs='+',
        metavar='INTERVAL', help='All intervals between two neighboring sites.')

    parser.add_argument('--mu',
        default=MU, type=int,
        metavar='MU', help='Expected shift in measurement relative to ref.')

    parser.add_argument('--sigma',
        default=SIGMA, type=int,
        metavar='SIGMA', help='Variance of shift in measurement relative to ref.')

    parser.add_argument('gap_name',
        metavar='GAP_NAME', help='The gaps name, give a name whatever you want.')

    parser.add_argument("-l", "--log", dest="log_level",
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='INFO',
        metavar='LOG_LEVEL', help="Set the logging level.")

    args=parser.parse_args()

    MU = args.mu
    SIGMA = args.sigma

    # with open(args.graph_file, 'rb') as fin:
        # nodes = pickle.load(fin)
        # sites = pickle.load(fin)

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
    logger.info('Start site ids (%d ids): %s', len(args.start_sites), args.start_sites)
    logger.info('End site ids (%d ids): %s', len(args.end_sites), args.end_sites)
    logger.info('Measurements (%d intervals): %s', len(args.intervals), args.intervals)
    

if __name__ == "__main__":
    main()