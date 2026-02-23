# Tetrahedra Untangler

This repository contains a 3D volumetric mesh untangler designed to fix inverted tetrahedra in 3D meshes. It is a direct reimplementation of the algorithm from [ssloy/invertible-maps](https://github.com/ssloy/invertible-maps), specifically adapted to be compatible with the **Geogram** geometry processing library.

## Overview

The core algorithm formulates mesh untangling as an optimization problem, adjusting vertex positions to maximize element quality and eliminate inverted elements (where the Jacobian determinant is $\le 0$). 

This implementation uses Geogram's built-in data structures (`GEO::Mesh`) and its HLBFGS optimizer to minimize the deformation energy.

## Key Differences & Features
While the core math remains similar to the original implementation, several modifications have been made:

* **Geogram Integration:** Fully relies on Geogram for mesh handling ( `GEO::Mesh`, `GEO::Optimizer`, ...).
* **Thread-Local OpenMP Gradients:** The global gradient calculation (`G`) is parallelized using thread-local storage (`G_local_store_`) to prevent race conditions without relying on expensive atomic operations.
* **Understandable Gradients:** The gradient calculations have been slightly refactored for clarity. If you want to understand the underlying math and derivation for the foldover-free map, refer to [this computational geometry write-up](https://baykeremm.github.io/computational-geometry/foldover/).

## Performance & OpenMP Tuning (Important)

Careful attention has been paid to how OpenMP is integrated. **More threads do not automatically mean faster execution**, especially on modern hybrid CPU architectures (e.g., Intel processors with Performance and Efficiency cores).

When calculating scaling speedup, you may notice diminishing returns or even performance degradation if threads are scheduled on E-cores. For example, on a 16-core machine with 6 P-cores, scaling efficiency drops off after 5-6 threads.

**Optimization Tip:** To get the best performance, restrict the execution of the process strictly to your Performance cores using `taskset`. 

```bash
# Example: If your P-cores are mapped to logical CPUs 0 through 5
taskset -c 0-5 ./your_executable
```
## Dependencies 

- Geogram 
- OpenMP

## Build and run 

To build the example: 
```bash
mkdir build 
cd build 
cmake .. 
make run
```

This will run the example included in the directory, the `domain_02.msh` is the original domain, then there is the same domain but deformed a bit to create 430 negative tetrahedra. 
You can inspect the result with `vorpaview out.geogram`.
