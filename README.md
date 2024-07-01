# Path Survival Probabilities using Zero-Suppressed Binary Decision Diagram

Novel reliability measures tailored for stochastic graphs with designated source vertices and failure-probability-weighted edges called Path Survival Probabilities that pertain to:

- A per-vertex reliability measure quantifying the average survival likelihood of single-source paths from a vertex to any source
- An overall reliability measure incorporating the graph density and the shortest distance to a source as regulating elements

## Usage

### Build

```
make
```

### PSZDD

```
./pszdd <edgelist> <sources> [ <options>... ]
```

#### Example

```
./pszdd Graphs/kobe.dat Graphs/kobe_s.dat
```

## File Format

### edgelist

Example:

```
1 2 0.9
1 3 0.8
2 3 0.5
```

`<edgelist>` specifies the edges of the network with their attributes. Each line starts with the two vertices corresponding to the endpoints of the edge. This is followed by the probability of non-failure for the edge. The edge length follow if the flag `len` is used.

### sources

Example:

```
1 3
```

`<sources>` specifies the source vertices in the network.

## Preprint

- [Path Survival Probabilities as Measures of Reliability for Lifeline Utility Networks](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4538362)

## Related Repositories

- [TdZdd Library](https://github.com/kunisura/TdZdd/)
- [PIZDD Repository](https://github.com/renzopereztan/PIZDD)
