# simpliciality_test

explore possible algorithms for simpliciality test

## Proof of concept
First things first, install the required libraries:
```python
pip install -r requirements.txt
```

To warm up, try:
```python
python is_simplicial.py -k datasets/00_degs.txt -s datasets/00_sizes.txt
```
By running the snippet several times, you may notice that the algorithm is not deterministic. 
Indeed, we are building the state space tree with the backtracking algorithm, 
and therefore sampling the ensemble that satisfies the joint degree sequence
(defined by `00_degs.txt` & `00_sizes.txt`). In other words, the joint sequence is simplicial!

The backtracking algorithm is an exhaustive search; 
that is, all feasible solutions are considered and it will always find the optimal solution. 
This means that, unfortunately, our sampling algorithm only works for very small systems.  

We still lack a proper existence test to check whether a joint sequence is simplicial. 
Luckily, when the sequence is simplicial, 
we do find a greedy deterministic algorithm that picks up an simplicial instance! 
This can work in fairly large inputs.

For example,   
```python
python is_simplicial.py -k datasets/01_degs.txt -s datasets/01_sizes.txt --greedy
```
The dataset is the [crime network dataset](https://github.com/jg-you/scm/blob/master/datasets/crime_facet_list.txt)
from the [Phys. Rev. E paper](https://doi.org/10.1103/PhysRevE.96.032312), having 551 nodes and 194 facets. 

Moreover, we find that the greedy algorithm can go with two directions,
and sometimes only the "backward direction" work.

For example,
```python
python is_simplicial.py -k datasets/02_degs.txt -s datasets/02_sizes.txt --greedy --no-forward
``` 

Interesting? I think it is!

## MISC notes (to clean up later)
* The graphical [Erdős–Gallai theorem](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93Gallai_theorem). 