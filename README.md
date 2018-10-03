# cluster_finder
This is an original project I set up for a basic course on algorithms and data structures.

This piece of software is an heuristic method to find clusters in facebook-like undirected graphs running in O(n log^2(n)).

Provided that node degrees are uniformly bounded we can partition our graph in cliques in linear time using a trivial heuristic that chooses first low-degree nodes in order to possibly build large cliques.

Then we construct an auxiliary weighted graph where cliques in the partition above are new nodes and arc (A,B) has weight equals to arc density between cliques A and B for each A, B. Now we run Kruscal algorithm and cluster together a fixed fraction of total nodes (and get another log(n) in time complexity), construct another auxiliary graph like we've just done and iterate till we are satisfied about the modularity of clusters obtained.
