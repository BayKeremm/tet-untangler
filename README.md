# Tetrahedra Untangler

This repository contains a 3D volumetric mesh untangler designed to fix inverted tetrahedra in 3D meshes. It is a direct reimplementation of the algorithm from [ssloy/invertible-maps](https://github.com/ssloy/invertible-maps), specifically adapted to be compatible with the **Geogram** geometry processing library.

## Overview

The core algorithm formulates mesh untangling as an optimization problem, adjusting vertex positions to maximize element quality and eliminate inverted elements (where the Jacobian determinant is $\le 0$). 

This implementation uses Geogram's built-in data structures (`GEO::Mesh`) and its HLBFGS optimizer to minimize the deformation energy. Boundary vertices are automatically locked during the optimization process to preserve the outer shape of the domain.

## Key Differences & Features

While the core math remains faithful to the original implementation, several modifications have been made for readability, framework integration, and specific use cases:

* **Geogram Integration:** Fully relies on Geogram for mesh handling, bounding volume hierarchies (`GEO::MeshFacetsAABB`), and optimization (`GEO::Optimizer`).
* **Thread-Local OpenMP Gradients:** The global gradient calculation (`G`) is parallelized using thread-local storage (`G_local_store_`) to prevent race conditions without relying on expensive atomic operations.
* **Understandable Gradients:** The gradient calculations have been slightly refactored for clarity. If you want to understand the underlying math and derivation for the foldover-free map, refer to [this computational geometry write-up](https://baykeremm.github.io/computational-geometry/foldover/).

## Surface Distance Penalty (Experimental)

The biggest functional addition to this codebase is a custom penalty term designed to minimize the distance of specific "interface" vertices to a target surface. 

* **How it works:** When enabled, it calculates the squared distance from flagged interface vertices to the closest facet on a tracking mesh using an AABB tree, adding this to the total energy.
* **Usage:** This feature is currently in development for a specific application. It can be entirely ignored/disabled by setting the macro `TUNE=0` in the source code, which reverts the behavior to a standard untangler.

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
