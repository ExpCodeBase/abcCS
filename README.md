# README #

This is the codebase of the Space-Efficient abcCS. 

### Graph format

- `graph.meta` contains the number of edges and the number of nodes in each layer. See `data/example/graph.meta`.

- `graph.e` contains all the edges. See `data/example/graph.e`.

### Usage

- Dominant coreness computation via delta: `./skyind -skyIndDelta path_of_datasets`.

- Query via B-Index: `./skyind -skyIndLv1 path_of_datasets alpha beta vq`.

- Parallel construction of B-Index: `./skyind -skyIndPar path_of_datasets thread_num`.


