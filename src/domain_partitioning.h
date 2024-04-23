#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
struct PartitionOperators {
    Eigen::SparseMatrix<double> R;
    Eigen::SparseMatrix<double> E;
};

std::vector<PartitionOperators>
overlapping_decomposition(const Eigen::MatrixXi &T, int num_nodes, int num_partitions, int overlap_layers);

std::vector<int> metis_partition(const Eigen::MatrixXi &T, int num_nodes, int num_partitions);

std::vector<std::vector<int>> create_overlapping_partitions(const Eigen::MatrixXi &T, int num_nodes, int num_partitions,
                                                            const std::vector<int> &element_partitions,
                                                            int overlap_layers);