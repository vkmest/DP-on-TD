# Vittorio Vičević (12347507)
# check if the graph is 3-colorable


import matplotlib.pyplot as plt
from collections import defaultdict
import itertools


def is_valid_coloring(graph, coloring):
    for u, v in graph.edges():
        if coloring[u] == coloring[v]:
            return False
    return True


def dp_3_colorability(graph, tree_decomposition):
    colors = [0, 1, 2]
    
    bags = list(tree_decomposition)
    dp_table = defaultdict(dict)

    for i, bag in enumerate(bags):
        for assignment in itertools.product(colors, repeat=len(bag)):
            coloring = dict(zip(bag, assignment))
            
            if is_valid_coloring(graph.subgraph(bag), coloring):
                dp_table[i][assignment] = coloring

    def merge(parent, child):
        parent_bag = list(bags[parent])
        child_bag = list(bags[child])
        
        intersection = set(parent_bag) & set(child_bag)
        new_table = {}

        for p_assignment, p_coloring in dp_table[parent].items():
            for c_assignment, c_coloring in dp_table[child].items():
                if all(p_coloring[v] == c_coloring[v] for v in intersection):
                    merged_coloring = p_coloring.copy()
                    merged_coloring.update(c_coloring)
                    new_table[p_assignment] = merged_coloring

        dp_table[parent] = new_table

    for i in range(len(bags) - 1, 0, -1):
        merge(i - 1, i)

    root_assignments = dp_table[0]
    if root_assignments:
        print("The graph is 3-colorable.")
        
        coloring = next(iter(root_assignments.values()))
        
        print("One valid 3-coloring solution:", coloring)
        return coloring
    else:
        print("The graph is not 3-colorable.")
        return None
