import sys


class CMap:
    def __init__(self, uid, positions, length):
        self.id = uid
        self.positions = positions
        self.length = length
        fragments = []
        if positions:
            fragments.append(positions[0])
            for i in range(1, len(positions)):
                fragments.append(positions[i] - positions[i - 1])
            fragments.append(length - positions[-1])
        self.fragments = fragments

    def __str__(self):
        return str((self.id, self.length, self.num_site))

    @property
    def num_site(self):
        return len(self.positions)

def read_file(file_name):
    result = {}
    with open(file_name) as fin:
        current_id = None
        positions = []
        for line in filter(lambda x: x[0] != '#', fin):
            tokens = line.split('\t', 6)
            uid, position = tokens[0], float(tokens[5])
            if not uid == current_id:
                if positions:
                    # Add a new cmap.
                    result[current_id] = CMap(current_id, positions[:-1], length)
                # Update loop variable.
                current_id, positions, length = uid, [position], float(tokens[1])
            else:
                # Continue building old cmap.
                positions.append(position)
        if positions:
            # Add last cmap.
            result[current_id] = CMap(current_id, positions[:-1], length)
    return result

def _test_read_file():
    input_file = sys.argv[1]
    for k, v in read_file(input_file).items():
        print(k, v)
        print(v.positions)
        print(v.fragments)


if __name__ == '__main__':
    _test_read_file()

