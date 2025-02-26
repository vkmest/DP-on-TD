from itertools import product

def is_3_colorable(G, tree_decomposition):
    root = next(iter(tree_decomposition.nodes())) #arbitrary root

    tree_parent = {}  
    tree_children = {node: [] for node in tree_decomposition.nodes()}

    def dfs(node, parent):      # parent-child relationships (directed tree)
        tree_parent[node] = parent
        for neighbor in tree_decomposition.neighbors(node):
            if neighbor != parent:
                tree_children[node].append(neighbor)
                dfs(neighbor, node)
    dfs(root, None)

    dp = {}     # DP table storing valid assignments for each bag

    def compute_dp(node):
        bag = tree_decomposition.nodes[node]['bag']
        dp[node] = {}

        # generate all valid color assignments for the bag
        valid_assignments = []
        for colors in product(range(3), repeat=len(bag)):
            # filters invalid assignments where 2 connected nodes get the same color
            if all(colors[i] != colors[j] for i in range(len(bag)) for j in range(i+1, len(bag)) if G.has_edge(bag[i], bag[j])):
                valid_assignments.append(colors)

        # children nodes are first
        for child in tree_children[node]:
            compute_dp(child)

        # we want to check if the assignment can exted to children
        for assignment in valid_assignments:
            assignment_pointers = {}
            extendable = True

            for child in tree_children[node]:
                child_bag = tree_decomposition.nodes[child]['bag']
                common = set(bag) & set(child_bag)      # shared nodes between bag and child bag
                found = False

                for child_assignment in dp[child]:  
                    # we want to check if child's valid assignment matches the parent's assignment on shared nodes
                    if all(assignment[bag.index(v)] == child_assignment[child_bag.index(v)] for v in common):
                        found = True
                        assignment_pointers[child] = child_assignment
                        break
                
                if not found:
                    extendable = False
                    break
            
            if extendable:  # if at least one child assignment is valid, store it 
                dp[node][assignment] = assignment_pointers

    compute_dp(root)

    if not dp[root]:
        return False, None  # not 3-colorable

    # reconstruct coloring (one valid assignment from the root bag)
    root_assignment = next(iter(dp[root]))
    coloring = {}

    def reconstruct(node, assignment):
        bag = tree_decomposition.nodes[node]['bag']
        for i, v in enumerate(bag):
            if v not in coloring:
                coloring[v] = assignment[i]

        for child in tree_children[node]:
            # recusively assign colors to nodes
            reconstruct(child, dp[node][assignment][child])

    reconstruct(root, root_assignment)
    return True, coloring
