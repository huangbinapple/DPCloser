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

    def reverse(self, site_len):
        self.fragments.reverse()
        self.positions = [self.length - self.positions[i] - site_len + 2
            for i in reversed(range(len(self.positions)))]  # 2 due to 1-index nature of `positions`.
        for ele in self.fragments:
            self.positions.append(self.positions[-1] + ele)

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

def compare_cmaps(cmaps_, cmaps, align_option=True):
    transform_arrays = {}
    for cmap_id, cmap_ in cmaps_.items():
        if len(cmap_.positions) != len(cmaps[cmap_id].positions):
            transform_array = []
            transform_arrays[cmap_id] = transform_array
            index = 0
            positions = cmaps[cmap_id].positions
            for position_ in cmap_.positions:
                if position_ == positions[index]:
                    transform_array.append(index)
                else:
                    if align_option:
                        transform_array.append(index)
                    else:
                        transform_array.append(index + 1)
                    index += 1
                index += 1
            assert len(transform_array) == len(cmap_.positions)
        else:
            transform_arrays[cmap_id] = list(range(len(cmap_.positions)))
    return transform_arrays

def _test_read_file():
    input_file = sys.argv[1]
    for k, v in read_file(input_file).items():
        print(k, v, len(v.positions))
        print(v.positions)
        print(v.fragments)


if __name__ == '__main__':
    _test_read_file()

