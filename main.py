# Vittorio Vicevic

import networkx as nx
from graph import read_graph, draw_graph
from tree_decomposition import compute_tree_decomposition
from colorability import is_3_colorable

if __name__ == "__main__":
    file_path = "graphs/graph12.txt"  # change this to your test file (graph00-graph13)

    G = read_graph(file_path)   # read graph
    draw_graph(G, title="Starting Graph")

    treewidth, tree_decomposition = compute_tree_decomposition(G)  # compute tw and td
    
    print("Tree Decomposition:")
    for i, bag in enumerate(tree_decomposition.nodes()):
        print(f"Bag {i}: {tree_decomposition.nodes[bag]['bag']}")
    print(f"Computed Treewidth: {treewidth}")  

    is_colorable, coloring = is_3_colorable(G, tree_decomposition) # check 3-colorability

    if is_colorable:
        print("The graph is 3-colorable.")
        draw_graph(G, coloring=coloring, title="3-Colorable Graph")
    else:
        print("The graph is NOT 3-colorable.")

    # compute the chromatic number using NetworkX's greedy algorithm
    greedy_coloring = nx.coloring.greedy_color(G, strategy="largest_first")  
    chromatic_number = max(greedy_coloring.values()) + 1  # colors are starting from 0 so we add 1
    print(f"Approximate Chromatic Number (largest_first): {chromatic_number}")

    dsatur_coloring = nx.coloring.greedy_color(G, strategy="DSATUR")
    dsatur_chromatic = max(dsatur_coloring.values()) + 1
    print(f"Approximate Chromatic Number (DSATUR): {dsatur_chromatic}")
