#pragma once
#include "geogram/basic/numeric.h"
#include "geogram/mesh/mesh.h"
#include "geogram/mesh/mesh_AABB.h"
#include "geogram/numerics/optimizer.h"
#include "untangler_utils.hpp"

using namespace untangler;

class Untangler3D
{
  public:
    Untangler3D(GEO::Mesh &reference_mesh, GEO::Mesh &deformed_mesh);

    void apply_scaling();
    void evaluate_jacobian(const double *x);
    double evaluate_energy(const double *x);

    static void StaticFunc(GEO::index_t n, double *X, double &F, double *G)
    {
        current_instance_->Func(n, X, F, G);
    }

    static void static_newiter_callback(GEO::index_t it,
                                        const double *x,
                                        double f,
                                        const double *g,
                                        double gnorm)

    //     Logging (current energy, norm of gradient, step size)
    //     Checking custom stopping criteria
    //     Updating visualizations or debug output.
    {
        if (current_instance_)
        {
        }

    }

    void Func(GEO::index_t n, double *x, double &f, double *g);
    UntangleResult go();
    void undo_scaling();
    void lock_boundary_vertices();
    void bbox(GEO::vec3 &min, GEO::vec3 &max);
    void initialize_ref_tets();

  private:
    GEO::Optimizer_var optimizer;
    int voted_orientation;
    GEO::Mesh &reference_mesh;
    GEO::Mesh &deformed_mesh;
    double theta = 0.5; // 0.7 it was
    // the energy is (1-theta)*(shape energy) + theta*(area energy)

    int maxiter = 10000;   // max number of outer iterations
    int bfgs_maxiter = 10; // max number of inner iterations 20 it was
    double bfgs_threshold = .01;
    static Untangler3D *current_instance_;

    int debug = 1; // verbose level

    std::vector<std::array<GEO::vec3, 4>> ref_tets;

    std::vector<double> X;        // current geometry, vertices copied here
    std::vector<bool> locked_;    // currently lock = boundary vertices

    std::vector<std::vector<double>> G_local_store_;

    std::vector<GEO::mat3> J; // per-tet Jacobian matrix
    std::vector<GEO::mat3> K; // per-tet dual basis: det J = dot J[i] * K[i]
    std::vector<double> det;  // per-tet determinant of the Jacobian matrix
    std::vector<GEO::vec3> original_positions;

    double eps; // regularization parameter, depends on min(jacobian)

    double detmin; // min(jacobian) over all tetrahedra
    int ninverted; // number of inverted tetrahedra

    GEO::vec3 bbmin_, bbmax_;
    const double boxsize = 1.5;
    const double shrink = 1.1;
};

