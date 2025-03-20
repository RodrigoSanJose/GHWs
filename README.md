# GHWs
A Sage package for computing the generalized Hamming weights (GHWs) and the relative generalized Hamming weights (RGHWs) of linear codes. For more details about the algorithms used and the implementation, please read the upcoming paper. 

Examples can be found in the help section of each function (function? will return this section). 

**How to use**: download GHWs.py and write `load("GHWs.py")` in Sage. Depending on your current directory, this may require specifying the path of the package, e.g., `load("/home/user/GHWs.py")`. 

**Recommendation**: many of the main functions in this package have a verbose argument. I recommend using verbose=True when doing heavy computations, since the output can be used to estimate whether the algorithm is going to finish in a reasonable time or not. Moreover, since this provides the current upper and lower bounds, if we know an upper or lower bound by any other method, this information may allow us to determine the GHWs before the algorithm finishes. 

**Small example**:
```python
load("GHWs.py")
C = codes.BinaryReedMullerCode(1,4)
hierarchy(C)
```
Output: [8, 12, 14, 15, 16]


**List of functions**: We provide several main functions, as well as some auxiliary functions for working with generalized Hamming weights. Some of the main functions may use a noticeable amount of RAM for larger finite field sizes. We provide a low memory version of these functions that is slower in general, but with a lower RAM usage. 
  - **Main functions**: 'GHW', 'hierarchy', 'RGHW', 'rhierarchy', 'higher_spectrum', 'rhigher_spectrum', 'matrix_supp', 'subspaces', 'num_subspaces', 'wei_duality'.
  - **Low memory main functions**: 'GHW_low_mem', 'hierarchy_low_mem', 'RGHW_low_mem', 'rhierarchy_low_mem', 'higher_spectrum_low_mem', 'rhigher_spectrum_low_mem'.
  - **Auxiliary functions**: 'vecwt', 'coltw', 'standard', 'is_cyclic', 'bch_bound', 'information'.

**Tests**: it is possible to test that the functions are working propertly by running the test test_GHWs.sage. This can be done by writing `sage test_GHWs.sage`. This requires to have GHWs.py in the same directory (or the corresponding path has to be specified in the first line of test_GHWs.sage). The test should take between 130s and 400s, depending on the processor. 

**Citation**: if you use this implementation for your research, please consider citing the upcoming paper.
