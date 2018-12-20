class Node:

    """Assembly graph node"""
    def __init__(self, uid, cov, seq):
        """
        Initialize a node.
        Arguments:
            uid(str): unique id of node. i.e. 34, 45, 34r, 'r' means this is a reversed version.
            cov(float): average coverage of base in this node.
            seq(str): the seqence.
        """
        self.uid = uid
        self.coverage = cov
        self.seq = seq
        self.children = []
    
    def __str__(self):
        return "{{{}, {}}}".format(self.uid, [node.uid for (node, overlap) in self.children])

    def __repr__(self):
        return self.__str__()

    @staticmethod
    def merge(node_list):
        current_node = node_list[0]
        result = [current_node.seq]
        for node in node_list[1:]:
            for child_node, overlap_len in current_node.children:
                if child_node == node:
                    result.append(node.seq[overlap_len:])
            current_node = node
        return ''.join(result)

    @property
    def length(self):
        return len(self.seq)

    def addChild(self, node, overlap):
        """Add a child to `self`.
        Arguments:
            node(Node): a child node.
        """
        self.children.append((node, overlap))

def _testNodeHash():
    node = Node('3', 34, 'AAAA')
    dict_ = {}
    dict_[node] = 3
    print(dict_)
    print('hash:', hash(node))

