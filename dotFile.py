class DotFile:

    def __init__(self, output_file):
        self.output_file = output_file

    def __enter__(self):
        self._f = open(self.output_file, 'w')
        self._f.write('digraph G {\n')
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._f.write('}\n')
        self._f.close()

    def write(self, string):
        self._f.write(string)

    def add_attributes(self, attribute_dict):
        self._f.write(', '.join(('{}="{}"'.format(key, value) for key, value in attribute_dict.items())))

    def new_line(self):
        self._f.write('\n')

    def add_node(self, node_name, attribute_dict=dict()):
        self._f.write('"{}"'.format(node_name))
        if attribute_dict:
            self.write(' [')
            self.add_attributes(attribute_dict)
            self.write(']')
        self.new_line()

    def add_edge(self, node_a, node_b, attribute_dict={}):
        self._f.write('"{}"->"{}"'.format(node_a, node_b))
        if attribute_dict:
            self.write(' [')
            self.add_attributes(attribute_dict)
            self.write(']')
        self.new_line()


def main():
    pass

def _test_dot_file():
    with DotFile('test.dot') as f:
        f.write('subgraph cluster1 {\n')
        nodes = ('a', 'b')
        for node in nodes:
            f.add_node(node)
        f.write('}\n')
        f.write('subgraph cluster2 {\n')
        nodes = ('d', 'c')
        for node in nodes:
            f.add_node(node)
        f.write('}\n')
        f.add_edge('a', 'd', {"dir": "back"})
        f.add_edge('b', 'c')
    

if __name__ == '__main__':
    _test_dot_file()
