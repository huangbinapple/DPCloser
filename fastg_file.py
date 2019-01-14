from Bio import SeqIO
from assembly_graph import Node
from spades_util import read_long_name, read_short_name

def build_assembly_graph(fastg_file, overlap):
    # A dict hold all nodes, indexed by their uid.
    nodes = {}
    # A dict hold tuple of child node id, indexed by uid.
    children_index = {}
    for seq_record in SeqIO.parse(fastg_file, "fasta"):
        # Add a new node.
        short_name, children_list = read_long_name(seq_record.id)
        uid, length, cov, is_reversed = read_short_name(short_name)
        if is_reversed:
            uid += 'r'
        assert length == len(seq_record)
        node = Node(uid, cov, str(seq_record.seq))
        nodes[uid] = node
        # Record childNode ids for this node.
        child_ids = []
        for short_name in children_list:
            child_id, _, _, is_reversed = read_short_name(short_name)
            if is_reversed:
                child_id += 'r'
            child_ids.append(child_id)
        children_index[uid] = child_ids
    # print(list(children_index.keys()))
    # Add chldren for node in `nodes`.
    for node_id, node in nodes.items():
        for child_id in children_index[node_id]:
            node.addChild(nodes[child_id], overlap)
    return nodes

def _test_build_assembly_graph():
    test_file = "datasets/assembly_graph.fastg"
    nodes = build_assembly_graph(test_file, 127)
    print(len(nodes))

if __name__ == '__main__':
    _test_build_assembly_graph()
