"""This file transform .camp file from bionano to .opt file used in soma and omacc."""
import getopt
import sys


def read_cmap_file(cmap_file):
    fin = open(cmap_file)
    site_positions = []
    for line in filter(lambda x: not x.startswith('#'), fin):
        tokens = line.split('\t')
        site_positions.append(float(tokens[5]))
    fin.close()
    return site_positions

def wirte_file(file_name, species_name, site_name, intervals, std, is_circular=False):
    # Transfer interval unit from base to kbase.
    intervals = map(lambda x: x/1000, intervals)
    # Transfer std unit from base to kbase.
    std /= 1000
    fout = open(file_name, 'w')
    # Write files.
    fout.write('\t'.join([species_name, 'contig00001', site_name,
        'CIRCULAR=' + ('true' if is_circular else 'false')]))
    fout.write('\n')
    for counter, interval in enumerate(intervals):
        fout.write('\t'.join([str(counter), str(interval), str(std)]))
        fout.write('\n')

def print_help():
    body = '<-s site name> <-n species name> [-d standard deviation] [-c] <cmap_file> [opt_file]'
    detail = """
    -c: cmap_file is circular."""
    print('python3 {} {}'.format(__file__, body))
    print(detail)

def main():
    is_circular = False
    std = 400
    interface = 's:n:d:ch'
    species_name, site_name = None, None
    options, args = getopt.getopt(sys.argv[1:], interface)
    for option, value in options:
        if option == '-h':
            print_help()
            sys.exit()
        elif option == '-s':
            site_name = value
        elif option == '-n':
            species_name = value
        elif option == '-d':
            std = int(value)
        elif option == '-c':
            is_circular = True
    if len(args) == 2:
        cmap_file, opt_file = args
    elif len(args) == 1:
        cmap_file = args[0]
        opt_file = '.'.join([cmap_file.split('.')[0], 'opt'])

    site_positions = read_cmap_file(cmap_file)
    # Turn positions into intervals.
    intervals = [site_positions[0]]
    for i in range(1, len(site_positions)):
        intervals.append(site_positions[i] - site_positions[i-1])
    # Write file.
    wirte_file(opt_file, species_name, site_name, intervals, std, is_circular=is_circular)


if __name__ == "__main__":
    main()
