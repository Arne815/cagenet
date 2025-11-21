
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <fstream>
#include <iostream>
#include <string>


bool load_off(const std::string& path, Eigen::Matrix<long double, Eigen::Dynamic, 3>& vertices, Eigen::Matrix<int, Eigen::Dynamic, 3>& faces)
{
    std::ifstream in(path);
    if (!in) return false;

    std::string header;
    in >> header;
    if (header != "OFF") return false;

    int nv, nf, ne;
    in >> nv >> nf >> ne;

    vertices.resize(nv, 3);
    faces.resize(nf, 3);

    for (int i = 0; i < nv; i++)
    {
        double x, y, z;
        in >> x >> y >> z;
        vertices(i, 0) = (long double)x;
        vertices(i, 1) = (long double)y;
        vertices(i, 2) = (long double)z;
    }

    for (int i = 0; i < nf; i++)
    {
        int n, a, b, c;
        in >> n >> a >> b >> c;
        faces(i, 0) = a;
        faces(i, 1) = b;
        faces(i, 2) = c;
    }

    return true;
}


void computeMVCForOneVertex(const Eigen::Matrix<long double, Eigen::Dynamic, 3>& cage_vertices,
                            const Eigen::Matrix<int, Eigen::Dynamic, 3>& cage_faces,
                            const Eigen::Matrix<long double, 1, 3>& eta,
                            Eigen::Matrix<long double, 1, Eigen::Dynamic>& weights
                            )
{
    long double epsilon = 0.00000001;
    long double sumWeights = 0.0;
    weights.setZero();

    unsigned int num_cage_vertices = cage_vertices.rows();
    unsigned int num_cage_faces = cage_faces.rows();

    Eigen::Vector<long double, Eigen::Dynamic> d(num_cage_vertices);
    d.setZero();
    const Eigen::Matrix<long double, Eigen::Dynamic, 3> u;

}


void computeMVC(const Eigen::Matrix<long double, Eigen::Dynamic, 3>& cage_vertices,
                const Eigen::Matrix<int, Eigen::Dynamic, 3>& cage_faces,
                const Eigen::Matrix<long double, Eigen::Dynamic, 3>& eta_m,
                Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>& phi)
{

    // std::cout << cage_vertices.rows() << "  " << cage_vertices.cols() << std::endl;
    // std::cout << eta_m.rows() << "  " << eta_m.cols() << std::endl;

    phi.resize(cage_vertices.rows(), eta_m.rows());

    // std::cout << phi.rows() << " " << phi.cols() << std::endl;

    Eigen::Matrix<long double, 1, Eigen::Dynamic> weights(cage_vertices.rows());

    for (int eta_vertex_idx = 0; eta_vertex_idx < eta_m.rows(); eta_vertex_idx++)
    {
        const Eigen::Matrix<long double, 1, 3> eta = eta_m.row(eta_vertex_idx);
        computeMVCForOneVertex(cage_vertices, cage_faces, eta_m.row(eta_vertex_idx), weights);
        phi.col(eta_vertex_idx) = weights;

        std::cout << phi(0,0) << "  "  << phi(0,1) << std::endl;
        std::cout << phi(1,0) << "  "  << phi(1,1) << std::endl;

        break;
    }

    return;
}


// test main for single calculation of MVC
int main() {

    std::string mesh_filepath = "../test_data/mesh000.off";
    std::string cage_filepath = "../test_data/mesh000_ec.off";
    
    std::string human_seg = "../test_data/scape_corrected.txt";

    Eigen::Matrix<long double, Eigen::Dynamic, 3> cage_vertices;
    Eigen::Matrix<int, Eigen::Dynamic, 3> cage_faces;

    if (!load_off(cage_filepath, cage_vertices, cage_faces))
        std::cout << "ERROR loading cage off file" << std::endl;

    Eigen::Matrix<long double, Eigen::Dynamic, 3> mesh_vertices;
    Eigen::Matrix<int, Eigen::Dynamic, 3> mesh_faces;

    if (!load_off(mesh_filepath, mesh_vertices, mesh_faces))
        std::cout << "ERROR loading mesh off file" << std::endl;

    Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> phi;
    
    computeMVC(cage_vertices, cage_faces, mesh_vertices, phi);

    return 0;
}
