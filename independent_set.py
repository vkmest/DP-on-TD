# Vittorio Vičević (12347507)
# max weighted independent set


from collections import defaultdict
import itertools


def is_independent_set(graph, subset):
    for u in subset:
        for v in subset:
            if u != v and graph.has_edge(u, v):
                return False
    return True


def dp_max_weighted_independent_set(graph, tree_decomposition, weights):
    bags = list(tree_decomposition)
    
    dp_table = defaultdict(list)

    for i, bag in enumerate(bags):
        bag_list = list(bag)
        
        for subset in itertools.chain.from_iterable(itertools.combinations(bag_list, r) for r in range(len(bag_list) + 1)):
            subset_set = set(subset)
            if is_independent_set(graph.subgraph(bag_list), subset_set):
                weight = sum(weights[node - 1] for node in subset_set)
                dp_table[i].append((weight, subset_set))

    def merge(parent, child):
        parent_bag = list(bags[parent])
        child_bag = list(bags[child])
        
        intersection = set(parent_bag) & set(child_bag)
        new_table = []

        for p_weight, p_set in dp_table[parent]:
            for c_weight, c_set in dp_table[child]:
                if p_set & intersection == c_set & intersection:
                    merged_set = p_set | c_set
                    
                    if is_independent_set(graph, merged_set):
                        merged_weight = sum(weights[node - 1] for node in merged_set)
                        new_table.append((merged_weight, merged_set))

        if new_table:
            dp_table[parent] = new_table

    for i in range(len(bags) - 1, 0, -1):
        merge(i - 1, i)

    root_result = max(dp_table[0], key=lambda x: x[0])
    
    print(f"Maximum Weighted Independent Set Weight: {root_result[0]}")
    print(f"Maximum Weighted Independent Set: {root_result[1]}")
    
    return root_result[1]
