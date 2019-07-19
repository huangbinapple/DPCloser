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
    parser = argparse.ArgumentParser()

    parser.add_argument('graph_file', metavar='INPUT_FILE', help='A pickle file.')
    parser.add_argument('--start_sites', required=True, nargs='+',
        metavar='START_SITE_ID', help='All possible start site ids.')
    parser.add_argument('--end_sites', required=True, nargs='+',
        metavar='END_SITE_ID', help='All possible end site ids.')
    parser.add_argument('--intervals', required=True, type=float, nargs='+',
        metavar='INTERVAL', help='All intervals between two neighboring sites.')
    parser.add_argument('--mu', default=0, type=int,
        metavar='MU', help='Expected shift in measurement relative to ref.')
    parser.add_argument('--sigma', default=300, type=int,
        metavar='SIGMA', help='Variance of shift in measurement relative to ref.')
    parser.add_argument('gap_name',
        metavar='GAP_NAME', help='The gaps name, give a name whatever you want.')
    parser.add_argument("-l", "--log", dest="log_level",
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='INFO',
        metavar='LOG_LEVEL', help="Set the logging level.")

    args=parser.parse_args()

    global MU; MU = args.mu
    global SIGMA; SIGMA = args.sigma

    # with open(args.graph_file, 'rb') as fin:
        # nodes = pickle.load(fin)
        # sites = pickle.load(fin)

    logging.basicConfig(filename=args.gap_name + '.log', filemode='w',
        level=getattr(logging, args.log_level))
    
    logging.info('Input command: %s', ' '.join(sys.argv))
    logging.info('Gap name: %s', args.gap_name)
    logging.info('SIGMA: %d', SIGMA)
    logging.info('Start site ids: %s', args.start_sites)
    logging.info('End site ids: %s', args.end_sites)
    logging.info('Measurements: %s', args.intervals)
    

if __name__ == "__main__":
    main()