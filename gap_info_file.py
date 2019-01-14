# Absolete.
import sys
import xmap_file


class Gap:
    def __init__(self, start_node_id, end_node_id, 
        start_site_index, start_site_position,
        end_site_index, end_site_position, intervals):
        self.start_node_id = start_node_id
        self.end_node_id = end_node_id
        self.start_site_index = start_site_index
        self.start_site_position = start_site_position
        self.end_site_index = end_site_index
        self.end_site_position = end_site_position
        self.intervals = intervals

    def write_to_file(self, fout):
        fout.write(' '.join((self.start_node_id, self.end_node_id)))
        fout.write('\n')
        fout.write(str(self.start_site_index))
        fout.write(', ')
        fout.write(str(self.start_site_position))
        fout.write('; ')
        fout.write(str(self.end_site_index))
        fout.write(', ')
        fout.write(str(self.end_site_position))
        fout.write('\n')
        fout.write(' '.join(map(str, self.intervals)))
        fout.write('\n')

    @staticmethod
    def read_from_lines(lines):
        for i in range(3):
            line = lines.pop().rstrip()
            if i == 0:
                start_node_id, end_node_id = line.split(' ')
            elif i == 1:
                tokens = line.split('; ')
                start_site_index, start_site_position = \
                    map(int, tokens[0].split(', '))
                end_site_index, end_site_position = \
                    map(int, tokens[1].split(', '))
            else:
                intervals = map(int, line.split(' '))
        gap = Gap(start_node_id, end_node_id, start_site_index,
            start_site_position, end_site_index, end_site_position,
            intervals)
        return gap

def read_file(file_name):
    gaps = []
    with open(file_name) as fin:
        line = fin.readline()
        lines = []
        while line:
            if not line.startswith('#'):
                lines.append(line)
                if len(lines) == 3:
                    gap = Gap.read_from_lines(lines)
                    gaps.append(gap)
                    lines.clear()
            line = fin.readline()
    return gaps

def write_file(gaps, file_name, comments):
    with open(file_name, 'w') as fout:
        for c in comments:
            fout.write('# ' + c + '\n')
        for g in gaps:
            g.write_to_file(fout)

def _test_write_file():
    comments = ['hello', 'This is a file for test.']
    gaps = []
    gap = Gap('n1', 'n2', 45, 2324, 1, 232, [434, 234])
    gaps.append(gap)
    gap = Gap('n3', 'n4', 4, 2435, 2, 2442, [4, 234])
    gaps.append(gap)
    write_file(gaps, 'test.gap_info', comments)

def _test_read_file():
    gaps = read_file('test.gap_info')
    print(len(gaps))
    for gap in gaps:
        print(gap.start_node_id, gap.end_node_id)

def print_help():
    body = '<xmap file> <gap info file>'
    print('python3 {} {}'.format(__file__, body))


def main():
    pass


if __name__ == '__main__':
    # _test_write_file()
    _test_read_file()
