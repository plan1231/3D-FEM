#pragma once
#include <Eigen/Dense>
#include <string>

// Function to efficiently read a Tetrahedral mesh from a .msh file (version 4 format) using modern C++
// Parameters:
// - filePath: Path to the .msh file
// - TV: Eigen::MatrixXd to store vertex coordinates
// - TT: Eigen::MatrixXi to store tetrahedra connectivity
// - SF: Eigen::MatrixXi to store surface triangles
// Returns:
// - bool: True if the mesh is successfully read, false otherwise
bool readTetMesh_msh4(const std::string &filePath,
                      Eigen::MatrixXd &TV, Eigen::MatrixXi &TT,
                      Eigen::MatrixXi &SF);