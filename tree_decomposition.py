import networkx as nx
from networkx.algorithms.approximation import treewidth_min_fill_in

def compute_tree_decomposition(G):
    treewidth, tree_decomposition = treewidth_min_fill_in(G)

    for bag in tree_decomposition.nodes():
        tree_decomposition.nodes[bag]['bag'] = sorted(list(bag))

    return treewidth, tree_decomposition
