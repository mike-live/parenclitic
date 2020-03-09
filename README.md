# Parenclitic Network Generalized Algorithm implementation

Parenclitic is a Python package which can effectively produce network represenatation from numeric data.

- [More about parenclitic](#more-about-parenclitic)
- [Installation](#installation)
- [Getting Started](#getting-started)
- [The Team](#the-team)
- [Acknowledgements](#acknowledgements)

## More About parenclitic

The main idea is consider pairwise feature planes and decide is there a connection between 2 features based on control and deviated groups.
So, we consider 2 groups: control and deviated. Group of deviated samples somehow differ from control samples.
And we interested in features which can identify those distinction.
Here 2 cases arises: subject can distinct by one feature or they can be separated only by 2 features rather then 1.
First, we identify and exclude features that can distinguish samples only by linear case.
Second, we identify pairs of features and construct graph representation of those pairwise connections.
One node of network is a feature, and edge characterizes deviation of subject from control group by those 2 features.

Next step is a metric computation of graphs and understanding of underlying network complexity. 
Those metrics can be used as reduction of dimensionality for further ML algorithms.

To deal with those things we develop parenclitic library.

Our package provides 3 main features:
1. Build, save and load parenclitic network.
2. Choose or create kernel to identify edges.
3. Compute network metrics based on python-igraph package.

## Installation

Parenclitic is available on PyPI. You can install it through pip:
```
pip install parenclitic
```
Dependencies:
1. NumPy
2. python-igraph
3. Pandas
4. sklearn
5. scipy

Please, carefully check that python-igraph is correctly installed.

## Getting started

First load data. We generate it for example.
```python
    import numpy as np
    num_samples = 100
    num_features = 30
    shift = 2
    X = np.random.randn(num_samples, num_features)
    y = np.random.randint(2, size = (num_samples, ))
    mask = np.array(y, np.int32)
    X[mask == 0, :] += shift
    mask[y == 0] = -1
```

X - data values with 100 samples each with 30 features. 
y - vector with features labels (0, 1) (int type)
mask - vector with -1 means control group, +1 means devated group, +2 means test group (int type)

For example we shifts data for control group twice of standard deviation and we expect almost complete networks.

There are some steps to run parenclitic

0. Import parenclitic library
```python
    import parenclitic
```

1. Make kernel that decides is there is link between those pairs for particular subject.
For example it is a PDF kernel with automatically defined threshold.
```python
    kernel = parenclitic.pdf_kernel()
```

2. On some datasets groups can be easily separated by only one feature. To exclude such features IG_filter can be applied.
```python
    pair_filter = parenclitic.IG_filter()
```
These excluding can help to distinguish pair-based deviation from one-feature deviation.

3. Make parenclitic model which uses chosen kernel and filter.
```python
    clf = parenclitic.parenclitic(kernel = kernel, pair_filter = pair_filter)
```

4. Fit data using 2 workers and number of feature pairs per worker is 1000.
```python
    clf.fit(X, y, mask, num_workers = 2, chunk_size = 1000)
```

5. Save graphs as tsv (tab-separated values). Or you can choose 'npz' as NumPy zipped file.
```python
    clf.save_graphs(gtype = 'csv')
```

Full example you can see in src/parenclitic_sample.ipynb

## Parallel computation

Parallel computation based on multiprocessing library and it can paralellize feature pairs over multiple processes.

## The Team

Parenclitic project started by Krivonosov Mikhail in 2018 in Lobachevsky State University based on many works of M. Zanin, A. Zaikin.

## Acknowledgements

This work was supported by the megagrant "Digital personalized medicine for healthy aging (CPM-aging): network analysis of Large multi-omics data to search for new diagnostic, predictive and therapeutic goals" ¹ 074-02-2018-330 (1).

