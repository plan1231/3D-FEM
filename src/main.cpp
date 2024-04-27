#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/readMESH.h>
#include <igl/readOBJ.h>
#include <igl/boundary_facets.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <cmath>
#include <vector>
#include <igl/colormap.h>
#include <igl/stb/write_image.h>
#include <igl/embree/EmbreeRenderer.h>
#include <filesystem>
#include "Timer.h"
#include "domain_partitioning.h"
#include <omp.h>
#include <Eigen/src/SparseLU/SparseLU.h>
#include <CLI/CLI.hpp>

struct CLIArgs
{
    int num_partitions; // default is based on number of processors
    int num_overlap_layers = 3;
    int num_steps = 50;
    double timestep = 0.1;
    double alpha = 5.0;
    bool usePreconditioner = false;
    bool renderImages = false;
    CLIArgs(int argc, char *argv[])
    {
        CLI::App app{"3D FEM Solver with (optional) ASM Preconditioner"};

        app.add_option("--overlap", num_overlap_layers, "number of layers of overlap between partitions in the decomposition")
            ->default_val(num_overlap_layers);

        app.add_option("--partitions", num_partitions, "number of partitions to decompose the domain into")
            ->default_val(omp_get_num_procs());

        app.add_option("--steps", num_steps, "number of time steps to simulate")
            ->default_val(num_steps);

        app.add_option("--timestep", timestep, "time step size")
            ->default_val(timestep);

        app.add_option("--alpha", alpha, "material alpha")
            ->default_val(alpha);

        app.add_flag("-p,--preconditioner", usePreconditioner,
                     "use ASM preconditioner for the solver")
            ->default_val(usePreconditioner);

        app.add_flag("-r,--render", renderImages,
                     "render images of the simulation at each time step")
            ->default_val(renderImages);

        try
        {
            app.parse(argc, argv);
        }
        catch (const CLI::ParseError &e)
        {
            exit(app.exit(e));
        }
    }
};

class ASMPreconditioner
{
public:
    using MatrixType = Eigen::SparseMatrix<double>;

    Eigen::VectorXd solve(const Eigen::VectorXd &rhs) const
    {
        return preconditioner * rhs;
    }

    Eigen::ComputationInfo info() const
    {
        return preconditioner_computed ? Eigen::Success : Eigen::NumericalIssue;
    }

    void loadPartitions(std::vector<PartitionOperators> &&operators)
    {
        this->operators = std::move(operators);
    }

    void compute(const MatrixType &A)
    {
        preconditioner.resize(A.rows(), A.cols());
        preconditioner.setZero();
#pragma omp parallel for
        for (int i = 0; i < operators.size(); ++i)
        {
            MatrixType localA = operators[i].R * A * operators[i].E;
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
            solver.compute(localA);
            MatrixType identity(localA.rows(), localA.cols());
            identity.setIdentity();
            MatrixType localInverse = solver.solve(identity);
#pragma omp critical
            {
                std::cout << "Computed preconditioner for partition " << i << std::endl;
                preconditioner += operators[i].E * localInverse * operators[i].R;
            }
        }
        preconditioner.makeCompressed();
        std::cout << "finished computing preconditioner" << std::endl;
        std::cout << (double)preconditioner.nonZeros() / preconditioner.rows() * preconditioner.cols() << " percent of the preconditioner is non-zero" << std::endl;
        
        preconditioner_computed = true;
    }

private:
    std::vector<PartitionOperators> operators; // restriction and extension operators for each partition
    bool preconditioner_computed = false;
    Eigen::SparseMatrix<double> preconditioner;
};

namespace fs = std::filesystem;

void tetrahedron_stiffness(Eigen::Matrix<double, 4, 4> &Ke, double volume,
                           const Eigen::Ref<const Eigen::RowVector4i> &element,
                           const Eigen::Ref<const Eigen::MatrixXd> V)
{

    // Extract vertices of the tetrahedron
    Eigen::Vector3d X0 = V.row(element(0)).transpose();
    Eigen::Vector3d X1 = V.row(element(1)).transpose();
    Eigen::Vector3d X2 = V.row(element(2)).transpose();
    Eigen::Vector3d X3 = V.row(element(3)).transpose();

    // Compute the matrix of differences
    Eigen::Matrix3d DeltaX;
    DeltaX.col(0) = X1 - X0;
    DeltaX.col(1) = X2 - X0;
    DeltaX.col(2) = X3 - X0;

    // Compute the inverse of DeltaX to get transformation matrix T
    Eigen::Matrix3d T = DeltaX.inverse();

    // Compute the gradients of the basis functions
    Eigen::Vector3d N0 = T * Eigen::Vector3d(-1, -1, -1);
    Eigen::Vector3d N1 = T * Eigen::Vector3d(1, 0, 0);
    Eigen::Vector3d N2 = T * Eigen::Vector3d(0, 1, 0);
    Eigen::Vector3d N3 = T * Eigen::Vector3d(0, 0, 1);

    // Set the local stiffness matrix Ke
    Ke(0, 0) = volume * N0.dot(N0);
    Ke(0, 1) = volume * N0.dot(N1);
    Ke(0, 2) = volume * N0.dot(N2);
    Ke(0, 3) = volume * N0.dot(N3);

    Ke(1, 0) = Ke(0, 1);
    Ke(1, 1) = volume * N1.dot(N1);
    Ke(1, 2) = volume * N1.dot(N2);
    Ke(1, 3) = volume * N1.dot(N3);

    Ke(2, 0) = Ke(0, 2);
    Ke(2, 1) = Ke(1, 2);
    Ke(2, 2) = volume * N2.dot(N2);
    Ke(2, 3) = volume * N2.dot(N3);

    Ke(3, 0) = Ke(0, 3);
    Ke(3, 1) = Ke(1, 3);
    Ke(3, 2) = Ke(2, 3);
    Ke(3, 3) = volume * N3.dot(N3);
}

void tetrahedron_mass(Eigen::Matrix<double, 4, 4> &Me, double volume)
{
    Me.setConstant(volume / 20.0);
    Me.diagonal().setConstant(volume / 10.0);
}

double tetrahedron_volume(const Eigen::Vector3d &X0, const Eigen::Vector3d &X1, const Eigen::Vector3d &X2, const Eigen::Vector3d &X3)
{
    Eigen::Vector3d AB = X1 - X0;
    Eigen::Vector3d AC = X2 - X0;
    Eigen::Vector3d AD = X3 - X0;
    double volume = std::fabs(AB.dot(AC.cross(AD))) / 6.0;
    return volume;
}

int main(int argc, char **argv)
{
    CLIArgs args = CLIArgs(argc, argv);

    Timer timer;
    timer.new_activity("calculate preconditioner");
    timer.new_activity("solve");
    timer.new_activity("build stiffness and mass matrices");

    // Load geometric data
    Eigen::MatrixXd surfaceV, V; // Vertices
    Eigen::MatrixXi T;           // Tetrahedrons (connectivity)
    Eigen::MatrixXi surfaceF, F; // Faces

    igl::readOBJ("../data/armadillo.obj", surfaceV, surfaceF);
    igl::copyleft::tetgen::tetrahedralize(surfaceV, surfaceF, "pq7.0Y", V, T, F);

    // igl::readMESH("../data/coarser_bunny.mesh", V, T, F);
    // igl::boundary_facets(T, surfaceF);
    // surfaceF = surfaceF.rowwise().reverse().eval();

    int Nv = V.rows(); // Number of vertices
    int Ne = T.rows(); // Number of elements

    std::cout << "building stiffness and mass matrices" << std::endl;
    Eigen::SparseMatrix<double> K(Nv, Nv);
    Eigen::SparseMatrix<double> M(Nv, Nv);
    timer.start(2);
#pragma omp parallel
    {
        Eigen::SparseMatrix<double> K_local(Nv, Nv);
        Eigen::SparseMatrix<double> M_local(Nv, Nv);
#pragma omp for
        for (int tet_idx = 0; tet_idx < Ne; ++tet_idx)
        {
            Eigen::Matrix<double, 4, 4> Ke, Me;

            double volume = tetrahedron_volume(V.row(T(tet_idx, 0)), V.row(T(tet_idx, 1)),
                                               V.row(T(tet_idx, 2)), V.row(T(tet_idx, 3)));

            auto indices = T.row(tet_idx);

            tetrahedron_stiffness(Ke, volume, indices, V);
            tetrahedron_mass(Me, volume);

            for (int ti = 0; ti < 4; ++ti)
            {

                for (int tj = 0; tj < 4; ++tj)
                {
                    K_local.coeffRef(indices(ti), indices(tj)) += Ke.coeff(ti, tj);
                    M_local.coeffRef(indices(ti), indices(tj)) += Me.coeff(ti, tj);
                }
            }
        }
#pragma omp critical
        {
            K += K_local;
            M += M_local;
        }
    }
    K.makeCompressed();
    M.makeCompressed();
    timer.stop();

    igl::embree::EmbreeRenderer er;
    er.set_mesh(V, surfaceF, true);
    Eigen::Matrix3d rot_matrix;
    // Specify rotation matrix:
    //     10 degrees around X axis
    //    180 degrees around Y axis
    //      4 degrees around Z axis
    rot_matrix = Eigen::AngleAxisd(10 * igl::PI / 180.0, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(igl::PI, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(4 * igl::PI / 180.0, Eigen::Vector3d::UnitZ());
    er.set_rot(rot_matrix);

    // create output directory if not exists
    if (args.renderImages)
    {
        fs::path dir = "result";
        if (!fs::exists(dir))
            fs::create_directory(dir);
    }

    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> buffer_R(1000, 1000);
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> buffer_G(1000, 1000);
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> buffer_B(1000, 1000);
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> buffer_A(1000, 1000);

    // Parameters
    Eigen::VectorXd c(Nv);

    for (int v = 0; v < Nv; v++)
    {
        c(v) = V.coeff(v, 0) > 0.1;
    }

    Eigen::SparseMatrix<double> A = M + args.timestep * args.alpha * K;

    std::cout << "Decomposing domain" << std::endl;
    auto operators = overlapping_decomposition(T, V.rows(), args.num_partitions, args.num_overlap_layers);

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, ASMPreconditioner> solver_asm;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;

    if (args.usePreconditioner)
    {
        std::cout << "Computing preconditioner with " << operators.size() << " subdomains" << std::endl;
        solver_asm.preconditioner().loadPartitions(std::move(operators));
        timer.start(0);
        solver_asm.compute(A);
        timer.stop();
    }

    solver.compute(A);

    std::cout << "Solving system" << std::endl;

    timer.print_on_exit(true);
    for (int step = 0; step < args.num_steps; ++step)
    {

        // Save result to a PNG
        if (args.renderImages)
        {
            er.set_data(c);
            er.render_buffer(buffer_R, buffer_G, buffer_B, buffer_A);
            std::ostringstream ss;
            ss << "result/" << step << ".png";
            igl::stb::write_image(ss.str(), buffer_R, buffer_G, buffer_B, buffer_A);
        }

        // solve the system for the next timestep
        Eigen::VectorXd b = M * c;
        timer.start(1);
        if (args.usePreconditioner)
        {
            c = solver_asm.solveWithGuess(b, c);
        }
        else
        {
            c = solver.solveWithGuess(b, c);
        }
        timer.stop();

        std::cout << step << std::endl
                  << "iterations: " << solver.iterations() << std::endl;
    }
    return 0;
}
