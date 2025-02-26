import networkx as nx
import matplotlib.pyplot as plt

def read_graph(file_path):
    G = nx.Graph()
    with open(file_path, 'r') as f:
        lines = f.readlines()
        n = int(lines[0].strip())  
        _ = list(map(int, lines[1].strip().split()))  # we are ignoring weights

        for line in lines[2:]:
            u, v = map(int, line.strip().split())
            G.add_edge(u, v)
    return G

def draw_graph(G, coloring=None, title="Graph"):
    pos = nx.spring_layout(G)  # layout for visualization -> it spaces the nodes naturally
    if coloring:
        color_map = {0: 'red', 1: 'green', 2: 'blue'}   # use 0, 1, 2 for colors
        node_colors = [color_map[coloring[node]] if node in coloring else 'gray' for node in G.nodes()]
    else:       # no coloring
        node_colors = 'lightblue'
    
    plt.figure()  # Create a new figure for each graph
    nx.draw(G, pos, with_labels=True, node_color=node_colors, edge_color='gray')
    plt.title(title)
    
    # Set the window title (Fixes "Figure 1" issue)
    plt.gcf().canvas.manager.set_window_title(title)

    plt.show()
