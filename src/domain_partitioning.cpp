#include "domain_partitioning.h"
#include <unordered_set>
#include <set>
#include <metis.h>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

std::vector<PartitionOperators>
overlapping_decomposition(const Eigen::MatrixXi &T, int num_nodes, int num_partitions, int overlap_layers) {

    std::vector<int> element_partitions = metis_partition(T, num_nodes, num_partitions);
    std::vector<std::vector<int>> overlapping_partitions = create_overlapping_partitions(T, num_nodes, num_partitions,
                                                                                       element_partitions, overlap_layers);
    
    std::vector<PartitionOperators> operators;

    // For each partition
    for (int p = 0; p < num_partitions; p++) {
        // restriction and extension operators


        // Use set instead of unordered_set to ensure consistent ordering. Not exactly necessary?
        std::set<int> nodes_set;
        for(int element: overlapping_partitions[p]){
            for(int j = 0; j < T.cols(); j++){
                nodes_set.insert(T(element, j));
            }
        }

        Eigen::SparseMatrix<double> R(nodes_set.size(), num_nodes);
        Eigen::SparseMatrix<double> E(num_nodes, nodes_set.size());
        // For each local node
        int i = 0;
        for (int global_node: nodes_set) {
            // Set the corresponding entry in the restriction operator to 1
            R.insert(i, global_node) = 1;
            i++;
        }

        E = R.transpose();

        operators.push_back(PartitionOperators{R, E});
    }

    return operators;
}

std::vector<int> metis_partition(const Eigen::MatrixXi &T, int num_nodes, int num_partitions){
    idx_t num_elements = T.rows();     // Number of elements
    idx_t num_corners = T.cols();      // Number of corners, assuming T is (elements x corners)

    // METIS requires arrays for the start of each element's connectivity list and the connectivity list itself
    std::vector<idx_t> eptr(num_elements + 1);
    std::vector<idx_t> eind(T.size());

    // Convert the Eigen connectivity matrix to METIS format
    for (idx_t i = 0; i < num_elements; i++) {
        eptr[i] = i * num_corners;
        for (idx_t j = 0; j < num_corners; j++) {
            // METIS expects 0-based indices
            eind[i * num_corners + j] = T(i, j);
        }
    }
    eptr[num_elements] = num_elements * num_corners;  // Last element of eptr is a sentinel

    // Variables needed for METIS
    idx_t edgecut;  // This will store the result of the edgecut
    std::vector<idx_t> epart(num_elements);  // Array to store the element partitioning
    std::vector<idx_t> npart(num_nodes);     // Array to store the node partitioning (may not be needed)

    // METIS options
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_CONTIG] = 1;  // Ensure the partitions are contiguous
    options[METIS_OPTION_NUMBERING] = 0; // C-style 0-based numbering

    idx_t ncommon = 3;
    // PartMeshDual partitions the mesh into k parts (elements)
    int ret = METIS_PartMeshDual(&num_elements, &num_nodes, eptr.data(), eind.data(), NULL, NULL, &ncommon, &num_partitions, NULL, options, &edgecut, epart.data(), npart.data());

    if (ret != METIS_OK) {
        throw std::runtime_error("METIS partitioning failed");
    }

    // Convert epart to int type before returning
    return {epart.begin(), epart.end()};
}

std::vector<std::vector<int>> create_overlapping_partitions(const Eigen::MatrixXi &T, int num_nodes, int num_partitions,
                                                            const std::vector<int> &element_partitions,
                                                            int overlap_layers) {

    // node -> element mapping
    std::vector<std::unordered_set<int>> node_to_element(num_nodes);

    for (int i = 0; i < T.rows(); ++i) {
        for (int j = 0; j < T.cols(); ++j) {
            node_to_element[T(i, j)].insert(i);
        }
    }

    // This will hold the final overlapping partitions
    std::vector<std::unordered_set<int>> overlapping_partitions(num_partitions);

    for (int i = 0; i < element_partitions.size(); ++i) {
        overlapping_partitions[element_partitions[i]].insert(i);
    }

    for (int layer = 0; layer < overlap_layers; ++layer) {
        std::vector<std::unordered_set<int>> new_elements(num_partitions); // Temporary sets to hold new elements to be added

        for (int p = 0; p < num_partitions; ++p) {
            for (int element: overlapping_partitions[p]) {
                // Look at each node in the element
                for (int j = 0; j < T.cols(); ++j) {
                    int node = T(element, j);
                    // Add all elements adjacent to this node to the new_elements set for this partition
                    for (int adj_element: node_to_element[node]) {
                        // Check if the element is not already in the partition to avoid duplicates
                        if (overlapping_partitions[p].find(adj_element) == overlapping_partitions[p].end()) {
                            new_elements[p].insert(adj_element);
                        }
                    }
                }
            }
        }

        for (int p = 0; p < num_partitions; ++p) {
            overlapping_partitions[p].insert(new_elements[p].begin(), new_elements[p].end());
        }
    }

    std::vector<std::vector<int>> overlapping_partitions_vec(num_partitions);
    for (int p = 0; p < num_partitions; ++p) {
        overlapping_partitions_vec[p].assign(overlapping_partitions[p].begin(), overlapping_partitions[p].end());
    }

    return overlapping_partitions_vec;
}
