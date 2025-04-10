# GHWs
A Sage package for computing the generalized Hamming weights (GHWs) and the relative generalized Hamming weights (RGHWs) of linear codes. For more details about the algorithms used and the implementation, please check the associated paper in https://doi.org/10.48550/arXiv.2503.17764. 

## How to use
Download GHWs.py and write `load("GHWs.py")` in Sage. Depending on your current directory, this may require specifying the path of the package, e.g., `load("/home/user/GHWs.py")`. 

## Recommendations
Many of the main functions in this package have a verbose argument. I recommend using `verbose=True` when doing heavy computations (for example, `hierarchy(C, verbose=True)`), since the output can be used to estimate whether the algorithm is going to finish in a reasonable amount of time or not. Moreover, since this provides the current upper and lower bounds, if upper or lower bounds are known by any other method, this information may allow to determine the GHWs before the algorithm finishes. 

## Small example
```python
load("GHWs.py")
C = codes.BinaryReedMullerCode(1,4)
hierarchy(C)
```
Output: [8, 12, 14, 15, 16]

## List of functions
We provide several main functions, as well as some auxiliary functions for working with generalized Hamming weights. Some of the main functions may use a noticeable amount of RAM for larger finite field sizes. We provide a low memory version of these functions that is slower in general, but with a lower RAM usage. 
  - **Main functions**: 'GHW', 'hierarchy', 'RGHW', 'rhierarchy', 'higher_spectrum', 'rhigher_spectrum', 'matrix_supp', 'subspaces', 'num_subspaces', 'wei_duality'.
  - **Low memory main functions**: 'GHW_low_mem', 'hierarchy_low_mem', 'RGHW_low_mem', 'rhierarchy_low_mem', 'higher_spectrum_low_mem', 'rhigher_spectrum_low_mem'.
  - **Auxiliary functions**: 'vecwt', 'colwt', 'standard', 'is_cyclic', 'bch_bound', 'information'.
    
Each function has a description text that can be accessed with function_name? (for example, `hierarchy?`). This description text explains what the function does, the parameters that it requires, the format of the output, and provides examples.

## Tests
It is possible to test that the functions are working propertly by running the test test_GHWs.sage (test_GHWs_low_mem.sage for the low memory functions). This can be done by writing `sage test_GHWs.sage`. This assumes the folder structure follows that of this repository. Otherwise, the line `load(path2)` from the test file has to be changed to specify the path of GHWs.py. The test should take between 130s and 500s, depending on whether the low memory functions are used or not and the processor's performance. The rest of the files are performance tests used to obtain the tables and graphs of the associated paper.

## Citation
If you use this implementation, please consider citing the paper https://doi.org/10.48550/arXiv.2503.17764 and/or this repository.

### Citation for the paper
```
@article{sanjoseGHWsPackage,
      journal={ArXiv 2503.17764},
      title={An algorithm for computing generalized {H}amming weights and the {S}age package {GHW}s}, 
      author={Rodrigo San-Jos\'{e}},
      year={2025},
      eprint={2503.17764},
      archivePrefix={arXiv},
      primaryClass={cs.IT}
}
```
### Citation for the repository
```
@Misc{githubGHWs,
    AUTHOR = {San-Jos\'{e}, Rodrigo},
    TITLE = {{GHWs}: A {S}age package for computing the generalized {H}amming weights of a linear code. {G}it{H}ub repository},
    year = "2025",
    howpublished = "Available online: \url{https://github.com/RodrigoSanJose/GHWs}"
}
```
