# 6_3_TREE

### read data ###
with open("../../inputs/stronghold/rosalind_tree.txt", "r") as f:
    n = int(f.readline().rstrip())
    adj_list = []
    for item in f.readlines():
        _append = tuple(map(int, item.rstrip().split()))
        adj_list.append(_append)

### calculation ###
## initial calculation
graph = [set(adj_list.pop(0))]
for adj in adj_list:
    for st_idx in range(len(graph)):
        st = graph[st_idx]
        include = False
        if (adj[0] in st) or (adj[1] in st):
            include = True
            graph[st_idx] = graph[st_idx].union(set(adj))
            break
        # if adj set not in graph, append to graph
        if st_idx == (len(graph) - 1):
            if include:
                continue
            else:  # not included
                graph.append(set(adj))
## calculate minimal edges for tree graph
_res = n
for i in graph:
    _res -= len(i)
res = (len(graph) + _res - 1)

### print ###
print(res)
