class DotFile:

    def __init__(self, outputFile):
        self.outputFile = outputFile

    def __enter__(self):
        self._f = open(self.outputFile, 'w')
        self._f.write('digraph G {\n')
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._f.write('}\n')
        self._f.close()

    def write(self, string):
        self._f.write(string)

    def addAttrs(self, attributeDict):
        self._f.write(', '.join(('{}="{}"'.format(key, value) for key, value in attributeDict.items())))

    def newLine(self):
        self._f.write('\n')

    def addNode(self, nodeName, attributeDict=dict()):
        self._f.write('"{}"'.format(nodeName))
        if attributeDict:
            self.write(' [')
            self.addAttrs(attributeDict)
            self.write(']')
        self.newLine()

    def addEdge(self, nodeA, nodeB, attributeDict={}):
        self._f.write('"{}"->"{}"'.format(nodeA, nodeB))
        if attributeDict:
            self.write(' [')
            self.addAttrs(attributeDict)
            self.write(']')
        self.newLine()


def main():
    pass

def _testDotFile():
    with DotFile('test.dot') as f:
        f.write('subgraph cluster1 {\n')
        nodes = ('a', 'b')
        for node in nodes:
            f.addNode(node)
        f.write('}\n')
        f.write('subgraph cluster2 {\n')
        nodes = ('d', 'c')
        for node in nodes:
            f.addNode(node)
        f.write('}\n')
        f.addEdge('a', 'd', {"dir": "back"})
        f.addEdge('b', 'c')
    

if __name__ == '__main__':
    _testDotFile()
