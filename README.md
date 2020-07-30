# simpliciality_test

explore possible algorithms for simpliciality test

## Proof of concept
First things first, install the required libraries:
```python
pip install -r requirements.txt
```

Now, try:
```python
python is_simplicial.py -k degs.txt -s sizes.txt
```

It should tell you that the joint degree sequence defined by `degs.txt` & `sizes.txt` is simplicial! 
However, currently, the code only works for very small systems.   