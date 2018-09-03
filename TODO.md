TODO


*Long term*
- refactor code
- adapt to run across multiple cores (Numba jit decorator with parallel=True flag)
- implement correlation length analysis
- ask users (ie, Richard!) to use github requests feature for better issue tracking
- write tests and use continuous integration (Travis)

*requests from users, ie: Richard*
- ULAM pattern: start in the middle and spiral outwards. With the option for starting at E/S/W/N and running clockwise from this (which would correspond to rotation at end of layer creation). Odd and even-squares will have different behaviour. Clarify defect positioning.
 - rectangle grid as an option

*Now implemented*
- Abbreviate sequences of multiple zeros : 7x0 2x1 6x0  DONE
- Input break points in same way as basic : 7x0 b 3x0 b DONE
- Create directory for a single run DONE
- Create input file DONE
- Overlay multiple layers (3rd / 4th) DONE
- Allow different patterns in each layer DONE
- Correct the issue with rotating - the break points are irrespective of the layer being done. The break points shouldn't rotate. DONE
- Output image with greater resolution  DONE 
