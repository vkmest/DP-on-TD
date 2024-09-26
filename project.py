# Vittorio Vičević (12347507)
# SD&A

import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms.approximation import treewidth_min_fill_in
import colorability
import independent_set


def read_graph_from_file(file_path):
    G = nx.Graph()
    with open(file_path, 'r') as f:
        lines = f.readlines()
        weights = list(map(int, lines[1].strip().split()))
        for line in lines[2:]:
            u, v = map(int, line.strip().split())
            G.add_edge(u, v)
    return G, weights


def draw_graph_colorability(G, coloring):
    pos = nx.spring_layout(G)
    color_map = {0: 'red', 1: 'green', 2: 'blue'}
    node_color = [color_map[coloring[node]] for node in G.nodes()]
    nx.draw(G, pos, with_labels=True, node_color=node_color, edge_color='gray')
    plt.title("3-Colorability of the Graph")
    plt.gcf().canvas.manager.set_window_title("3-Colorability of the Graph")
    plt.show()


def draw_graph_max_independent_set(G, independent_set):
    pos = nx.spring_layout(G)
    node_color = ['orange' if node in independent_set else 'lightblue' for node in G.nodes()]
    nx.draw(G, pos, with_labels=True, node_color=node_color, edge_color='gray')
    plt.title("Maximum Weighted Independent Set")
    plt.gcf().canvas.manager.set_window_title("Maximum Weighted Independent Set")
    plt.show()


def draw_graph_min_vertex_cover(G, vertex_cover):
    pos = nx.spring_layout(G)
    node_color = ['green' if node in vertex_cover else 'lightblue' for node in G.nodes()]
    nx.draw(G, pos, with_labels=True, node_color=node_color, edge_color='gray')
    plt.title("Minimum Vertex Cover")
    plt.gcf().canvas.manager.set_window_title("Minimum Vertex Cover")
    plt.show()


def draw_graph(G, vertex_set=None, coloring=None, title="Graph"):
    pos = nx.spring_layout(G)
    node_color = 'lightblue'
    if vertex_set:
        node_color = ['orange' if node in vertex_set else 'lightblue' for node in G.nodes()]
    elif coloring:
        color_map = {0: 'red', 1: 'green', 2: 'blue'}
        node_color = [color_map[coloring[node]] for node in G.nodes()]
    nx.draw(G, pos, with_labels=True, node_color=node_color, edge_color='gray')
    plt.title(title)
    plt.gcf().canvas.manager.set_window_title(title)
    plt.show()


if __name__ == "__main__":
    for i in range(11):
        graph_file_path = f'graphs/graph{i:02}.txt'
        print(f"\nProcessing {graph_file_path}\n")
        
        graph, weights = read_graph_from_file(graph_file_path)
        
        draw_graph(graph, title=f"Starting Graph {graph_file_path}")
        treewidth, tree_decomposition = treewidth_min_fill_in(graph)
        
        print("Tree Decomposition:")
        for i, bag in enumerate(tree_decomposition):
            print(f"Bag {i}: {sorted(list(bag))}")
        print(f"Computed treewidth: {treewidth}")
        
        coloring = colorability.dp_3_colorability(graph, tree_decomposition)
        if coloring:
            draw_graph_colorability(graph, coloring)
            plt.pause(0.1)
        
        independent_set_result = independent_set.dp_max_weighted_independent_set(graph, tree_decomposition, weights)
        draw_graph_max_independent_set(graph, independent_set_result)
        plt.pause(0.1)
        
        vertex_cover_result = set(graph.nodes()) - set(independent_set_result)
        draw_graph_min_vertex_cover(graph, vertex_cover_result)
