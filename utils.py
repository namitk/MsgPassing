from union_find import union_find 

def min_spanning_tree(G):
    """ Minimum spanning tree using UF data structure and Kruskal algorithm """
    edges = []
    uf = union_find()
    for i in range(len(G)):
        for j in range(len(G[i])):
            edges.append((G[i][j],i,j))
    tree = []
    C={}
    R={}
    for u in range(len(G)):
        C[u]=u
    for u in range(len(G)):
        R[u]=0
    for W,u,v in sorted(edges):
        if uf.find(C, u) != uf.find(C, v):
            tree.append((u, v)) #, set(self.max_cliques[u]).intersection(set(self.max_cliques[v]))))
            uf.union(C, R, u, v)
    return tree

def DFS(G, root):
    visited = [False for x in G]
    dfs = [{'send':[], 'recv':[]} for x in G]
    stack = [root]
    while stack:
        cur = stack.pop()
        if not visited[cur]:
            for x in G[cur]:
                if not visited[x]:
                    dfs[x]['send'].append(cur)
                    dfs[cur]['recv'].append(x)
                stack.append(x)
            visited[cur]=True
    return dfs
