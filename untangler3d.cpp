#include "untangler3d.hpp"
#include "geogram/basic/geometry.h"
#include "geogram/basic/numeric.h"
#include "untangler_utils.hpp"
#include <algorithm>
#include <cmath>
#include <geogram/mesh/mesh_AABB.h>

#define EPS_FROM_THE_THEOREM 0

using namespace GEO;

using namespace untangler;

Untangler3D *Untangler3D::current_instance_ = nullptr;

Untangler3D::Untangler3D(Mesh &reference_mesh, Mesh &deformed_mesh)

    : reference_mesh(reference_mesh), deformed_mesh(deformed_mesh),
      ref_tets(reference_mesh.cells.nb()), X(reference_mesh.vertices.nb() * 3),
      locked_(reference_mesh.vertices.nb(), false),
      J(reference_mesh.cells.nb()), K(reference_mesh.cells.nb()),
      det(reference_mesh.cells.nb()),
      original_positions(reference_mesh.vertices.nb())

{
    for (index_t v = 0; v < original_positions.size(); v++)
    {
        original_positions[v] = reference_mesh.vertices.point(v);
    }

    voted_orientation = 1;
    apply_scaling(); // populates X scales to_track

    lock_boundary_vertices();

    initialize_ref_tets();

    int max_threads = omp_get_max_threads();
    G_local_store_.resize(max_threads);
    for (auto &vec : G_local_store_)
    {
        vec.resize(reference_mesh.vertices.nb() * 3, 0.0);
    }

    optimizer = GEO::Optimizer::create("HLBFGS");
    optimizer->set_verbose(true);
    optimizer->set_funcgrad_callback(StaticFunc);
    optimizer->set_newiteration_callback(static_newiter_callback);
    // Sets the gradient‑norm tolerance for convergence: when ∥∇f(X)∥ drops
    // below bfgs_threshold, the optimizer considers the solution converged and
    // stops.
    optimizer->set_epsg(bfgs_threshold);
    // Sets the maximum iteration count. The optimizer will not perform more
    // than bfgs_maxiter iterations, even if the gradient tolerance is not yet
    // met.
    optimizer->set_max_iter(bfgs_maxiter);
    // Tells the optimizer how many variables it must optimize
    optimizer->set_N(X.size());
    // Sets the memory parameter MM of the L‑BFGS part of HLBFGS. The optimizer
    // will store the last M pairs of step and gradient‑difference vectors to
    // approximate the inverse Hessian.
    optimizer->set_M(5);
}

void Untangler3D::apply_scaling()
{
    bbox(bbmin_, bbmax_);
    double maxside = std::max(
        {bbmax_.x - bbmin_.x, bbmax_.y - bbmin_.y, bbmax_.z - bbmin_.z});

    auto s = boxsize / (shrink * maxside);
    auto t = vec3(1, 1, 1) * boxsize / 2;
    auto c = (bbmax_ + bbmin_) / 2.;
    for (vec3 &p : deformed_mesh.vertices.points())
    {

        p = (p - c) * s + t;
    }

    for (int v : deformed_mesh.vertices)
    {
        for (int d : range(3))
        {
            X[v * 3 + d] = deformed_mesh.vertices.point(v).data()[d];
        }
    }
}

void Untangler3D::undo_scaling()
{
    double maxside = std::max(
        {bbmax_.x - bbmin_.x, bbmax_.y - bbmin_.y, bbmax_.z - bbmin_.z});

    auto si = shrink / boxsize * maxside;
    auto t = vec3(1, 1, 1) * boxsize / 2;
    auto c = (bbmax_ + bbmin_) / 2.;

    for (int v : deformed_mesh.vertices)
    {
        vec3 p = {X[v * 3 + 0], X[v * 3 + 1], X[v * 3 + 2]};
        deformed_mesh.vertices.point(v) = (p - t) * si + c;
    }
}

void Untangler3D::lock_boundary_vertices()
{
    for (auto c : deformed_mesh.cells)
    {
        for (int lf : range(4))
        {
            if (deformed_mesh.cells.adjacent(c, lf) == NO_CELL)
            {
                for (auto lv : range(3))
                {
                    auto v = deformed_mesh.cells.facet_vertex(c, lf, lv);
                    locked_[v] = true;
                }
            }
        }
    }
}

void Untangler3D::bbox(vec3 &min, vec3 &max)
{
    min = max = deformed_mesh.vertices.point(0);
    for (vec3 const &p : deformed_mesh.vertices.points())
    {
        for (int d = 0; d < 3; d++)
        {
            min[d] = std::min(min[d], p[d]);
            max[d] = std::max(max[d], p[d]);
        }
    }
}

void Untangler3D::initialize_ref_tets()
{
    static const int face_vert[4][3] = {
        {1, 3, 2}, {0, 2, 3}, {0, 3, 1}, {0, 1, 2}};

    // --- 1. Calculate Target Volume (Stiffness Normalization) ---
    // We calculate the average volume of valid cells using original positions
    double total_vol = 0.0;
    int valid_cells = 0;

    for (auto c : deformed_mesh.cells)
    {
        vec3 p_original[4];
        for (int lv = 0; lv < 4; ++lv)
        {
            p_original[lv] =
                original_positions[deformed_mesh.cells.vertex(c, lv)];
        }
        double vol = std::abs(Geom::tetra_signed_volume(
            p_original[0], p_original[1], p_original[2], p_original[3]));
        if (vol > 1e-20)
        {
            total_vol += vol;
            valid_cells++;
        }
    }
    double target_ref_vol = (valid_cells > 0) ? (total_vol / valid_cells) : 1.0;

    // --- 2. Determine Orientation ---
    int n_neg = 0;
    for (auto c : deformed_mesh.cells)
    {
        vec3 p[4];
        for (int lv = 0; lv < 4; ++lv)
            p[lv] =
                deformed_mesh.vertices.point(deformed_mesh.cells.vertex(c, lv));
        if (Geom::tetra_signed_volume(p[0], p[1], p[2], p[3]) < 0)
            n_neg++;
    }
    std::cout << "Number of negative tets: " << n_neg << std::endl;
    voted_orientation = (n_neg > (int)deformed_mesh.cells.nb() / 2) ? -1 : 1;

    // --- 3. Build Reference Tetrahedra ---
    for (auto c : deformed_mesh.cells)
    {
        vec3 p[4];

        for (int lv = 0; lv < 4; ++lv)
        {
            index_t v = deformed_mesh.cells.vertex(c, lv);
            p[lv] = original_positions[v];
        }

        double current_vol =
            std::abs(Geom::tetra_signed_volume(p[0], p[1], p[2], p[3]));

        double s = std::cbrt(target_ref_vol / current_vol);

        // Apply scaling around centroid
        vec3 center = (p[0] + p[1] + p[2] + p[3]) / 4.0;
        for (int i = 0; i < 4; ++i)
        {
            p[i] = center + (p[i] - center) * s;
        }

        for (int lf : range(4))
        {
            int v0 = face_vert[lf][0];
            int v1 = face_vert[lf][1];
            int v2 = face_vert[lf][2];

            vec3 e0 = p[v1] - p[v0];
            vec3 e1 = p[v2] - p[v0];

            // Normalize by target_ref_vol so all tets have equal "strength"
            ref_tets[c][lf] = voted_orientation * (cross(e0, e1) / 2.0) /
                              (3.0 * target_ref_vol);
        }
    }
}

void Untangler3D::evaluate_jacobian(const double *X)
{
    detmin = std::numeric_limits<double>::max();
    ninverted = 0;
#pragma omp parallel for reduction(min : detmin) reduction(+ : ninverted)
    for (index_t c = 0; c < deformed_mesh.cells.nb(); c++)
    {
        // Compute J (Jacobian matrix)
        mat3 &J = this->J[c];
        J.load_zero();

        // Note: ref_tets[c][i] and mesh_.cells.vertex(c, i)
        // are ok, we are not mixing facet and vertex indexing
        // since ref_tets[c][i] is del lambda_i stored there
        // (facet across vertex i)
        // Remember: all 4 vertices contribute
        for (int i = 0; i < 4; i++)
        {
            int v = deformed_mesh.cells.vertex(c, i);
            double ui0 = X[3 * v + 0];
            double ui1 = X[3 * v + 1];
            double ui2 = X[3 * v + 2];

            for (int k = 0; k < 3; k++)
            {
                J(0, k) += ref_tets[c][i][k] * ui0;
                // J(0,k) += partial(lambda_i / xk) * ui0
                J(1, k) += ref_tets[c][i][k] * ui1;
                // J(1,k) += partial(lambda_i / xk) * ui1
                J(2, k) += ref_tets[c][i][k] * ui2;
                // J(2,k) += partial(lambda_i / xk) * ui2
            }
        }

        det[c] = GEO::det(J);
        detmin = std::min(detmin, det[c]);
        ninverted += (det[c] <= 0);
        // Cofactor (dual basis) matrix K
        // if you switch i,j below, you compute the adj(J) = Cof(J)^T
        // then you need to get the row of b_col below in Func
        // below computes the Cof(J) and not adj(J)
        mat3 &K = this->K[c];

        // column 0
        K(0, 0) = J(1, 1) * J(2, 2) - J(1, 2) * J(2, 1);
        K(1, 0) = J(0, 2) * J(2, 1) - J(0, 1) * J(2, 2);
        K(2, 0) = J(0, 1) * J(1, 2) - J(0, 2) * J(1, 1);

        // column 1
        K(0, 1) = J(1, 2) * J(2, 0) - J(1, 0) * J(2, 2);
        K(1, 1) = J(0, 0) * J(2, 2) - J(0, 2) * J(2, 0);
        K(2, 1) = J(0, 2) * J(1, 0) - J(0, 0) * J(1, 2);

        // column 2
        K(0, 2) = J(1, 0) * J(2, 1) - J(1, 1) * J(2, 0);
        K(1, 2) = J(0, 1) * J(2, 0) - J(0, 0) * J(2, 1);
        K(2, 2) = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
    }
}

double Untangler3D::evaluate_energy(const double *X)
{
    evaluate_jacobian(X);
    double E = 0;
#pragma omp parallel for reduction(+ : E)
    for (index_t c = 0; c < deformed_mesh.cells.nb(); c++)
    {
        double chi_ = chi(eps, det[c]);
        auto a = J[c];
        vec3 c1(a(0, 0), a(1, 0), a(2, 0));
        vec3 c2(a(0, 1), a(1, 1), a(2, 1));
        vec3 c3(a(0, 2), a(1, 2), a(2, 2));
        double f =
            (dot(c1, c1) + dot(c2, c2) + dot(c3, c3)) / pow(chi_, 2.0 / 3.0);

        double g = (1 + square(det[c])) / chi_;
        E += (1 - theta) * f + theta * g;
    }

    return E;
}

void Untangler3D::Func(index_t n, double *X, double &F, double *G)
{

    int nthreads = omp_get_max_threads();

// 1. Reset Global Gradient
// std::fill is sequential; parallelize it for large N
#pragma omp parallel for schedule(static)
    for (index_t i = 0; i < n; ++i)
        G[i] = 0.0;
    F = evaluate_energy(X);
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::fill(G_local_store_[tid].begin(), G_local_store_[tid].end(), 0.0);
        double *Gl = G_local_store_[tid].data();

#pragma omp for schedule(static)
        for (index_t t = 0; t < deformed_mesh.cells.nb(); t++)
        {
            mat3 &a = J[t]; // Jacobian J
            mat3 &b = K[t]; // Cofactor K

            double c1 = chi(eps, det[t]);
            double c2 = pow(c1, 2.0 / 3.0);
            double c3 = chi_deriv(eps, det[t]);

            // ---- compute f = tr JT J ----
            vec3 col1(a(0, 0), a(1, 0), a(2, 0));
            vec3 col2(a(0, 1), a(1, 1), a(2, 1));
            vec3 col3(a(0, 2), a(1, 2), a(2, 2));
            double f =
                (dot(col1, col1) + dot(col2, col2) + dot(col3, col3)) / c2;

            double g = (1.0 + det[t] * det[t]) / c1;

            // ---- Assemble Stress Tensor P ----
            mat3 P;

            // Loop over columns of J (reference coordinates j=0,1,2)
            // to compute columns of P
            for (int j = 0; j < 3; j++)
            {
                // Column j of J
                vec3 a_col(a(0, j), a(1, j), a(2, j));
                // Column j of K (Cofactor)
                vec3 b_col(b(0, j), b(1, j), b(2, j));

                // vec3 b_col(b(j, 0), b(j, 1), b(j, 2));
                // this is because if we computed the adj(J) and that is
                // Cofactor^T that is why we need to use the rows of K, this is
                // what I missed that cus i was computinmg adj(J) above in
                // ref_tets and now Cof(J) is computed

                // df/dJ_col_j (Part of stress tensor)
                vec3 dfda =
                    a_col * (2.0 / c2) - b_col * ((2.0 * f * c3) / (3.0 * c1));

                // dg/dJ_col_j
                vec3 dgda = b_col * ((2.0 * det[t] - g * c3) / c1);

                vec3 dEda = dfda * (1.0 - theta) + dgda * theta;

                // Store in matrix P (column j)
                P(0, j) = dEda.x;
                P(1, j) = dEda.y;
                P(2, j) = dEda.z;
            }

            // ---- Distribute Force to 4 vertices ----
            // Gradient formula: G_v = P * grad_lambda_v
            for (int i = 0; i < 4; i++)
            {
                index_t v = deformed_mesh.cells.vertex(t, i);
                if (locked_[v])
                    continue;

                // ref_tets[t][i] is grad_lambda (vector of size 3)
                vec3 grad_lambda = ref_tets[t][i];

                // Matrix-Vector multiplication: P * grad_lambda
                // P(row, col) * vec(col)
                double gx = P(0, 0) * grad_lambda.x + P(0, 1) * grad_lambda.y +
                            P(0, 2) * grad_lambda.z;
                double gy = P(1, 0) * grad_lambda.x + P(1, 1) * grad_lambda.y +
                            P(1, 2) * grad_lambda.z;
                double gz = P(2, 0) * grad_lambda.x + P(2, 1) * grad_lambda.y +
                            P(2, 2) * grad_lambda.z;

                // Accumulate

                Gl[v * 3 + 0] += gx;
                Gl[v * 3 + 1] += gy;
                Gl[v * 3 + 2] += gz;
            }
        }
    }

#pragma omp parallel for schedule(static)
    for (index_t i = 0; i < n; i++)
    {
        double sum = 0.0;
        for (int tid = 0; tid < nthreads; tid++)
            sum += G_local_store_[tid][i];
        G[i] = sum;
    }
}

UntangleResult Untangler3D::go()
{
    current_instance_ = this;
#pragma omp parallel
    {
#pragma omp single
        {
            std::cout << "OMP threads = " << omp_get_num_threads() << std::endl;
        }
    }
    evaluate_jacobian(X.data());

#if EPS_FROM_THE_THEOREM
    eps = 1.;
#else
    double e0 = 1e-3;
#endif

    //
    for (int iter = 0; iter < maxiter; iter++)
    {

#if !EPS_FROM_THE_THEOREM && !NO_EPS
        if (iter && iter % 10 == 0 && e0 > 1e-8)
            e0 /= 2.;
        eps = detmin > 0 ? e0 : std::sqrt(square(e0) + 0.04 * square(detmin));
#endif

        if (debug > 0)
        {
            std::cerr << "iteration #" << iter << std::endl;
            std::cerr << "E: " << evaluate_energy(X.data()) << " eps: " << eps
                      << " detmin: " << detmin << " ninv: " << ninverted
                      << std::endl;
        }

        double E_prev = evaluate_energy(X.data());
        // {
        //     TimerOMP t("optimize");
        optimizer->optimize(X.data());
        // }
        double E = evaluate_energy(X.data());

#if EPS_FROM_THE_THEOREM
        double sigma = std::max(1. - E / E_prev, 1e-1);
        if (detmin >= 0)
            eps *= (1 - sigma);
        else
            eps *= 1 - (sigma * std::sqrt(square(detmin) + square(eps))) /
                           (std::abs(detmin) +
                            std::sqrt(square(detmin) + square(eps)));

#endif

        if (detmin > 0 && std::abs(E_prev - E) / E < 1e-5)
        {
            break;
        }
    }

    undo_scaling();
    return UntangleResult(ninverted > 0, evaluate_energy(X.data()));
}
