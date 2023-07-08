# Path Survival Probabilities using Zero-Suppressed Binary Decision Diagram

Novel reliability measures tailored for graphs with designated source vertices and failure-probability-weighted edges called Path Survival Probabilities that pertain to:

- A per-vertex reliability measure quantifying the average survival likelihood of single-source paths from a vertex to any source
- An overall reliability measure incorporating the graph density and the shortest distance to a source as regulating elements

## Usage

### Build

```
make
```

### PSZDD

```
./pszdd <edgelist> <probabilities> <sources> [ <options>... ]
```

#### Example

```
./pszdd Graphs/kobe_graph.dat Graphs/kobe_p.dat Graphs/kobe_t.dat
```

## File Format

### edgelist

Example:

```
1 2
1 3
2 3
```

`<edgelist>` is an edge list of the network. Each line consists of two vertices corresponding to the endpoints of an edge in the network.

### probabilities

Example:

```
0.9 0.8 0.7
```

`<probabilities>` contains the probabilities of non-failure for each edge in the `<edgelist>`. The order of the probabilities corresponds to the order of the edge list.

### sources

Example:

```
1 3
```

`<sources>` specifies the source vertices in the network.

## Related Repositories

- [TdZdd Library](https://github.com/kunisura/TdZdd/)
- [PIZDD Repository](https://github.com/renzopereztan/PIZDD)
