#include "readTetMesh.h"
#include <igl/boundary_facets.h>
#include <fstream>

bool readTetMesh_msh4(const std::string &filePath,
                      Eigen::MatrixXd &TV, Eigen::MatrixXi &TT,
                      Eigen::MatrixXi &SF)
{
    std::ifstream fileStream(filePath); // Open file for reading
    if (!fileStream.is_open())
    {
        return false; // Unable to open file, return false
    }

    std::string line;
    while (std::getline(fileStream, line))
    { // Read file line by line
        // Check for section indicating node coordinates
        if (line.find("$Nodes") != std::string::npos)
        {
            int bypass, numNodes;
            fileStream >> bypass >> numNodes >> bypass >> bypass; // Total nodes
            TV.resize(numNodes, 3);                               // Resize the matrix to accommodate vertices
            assert(numNodes > 0);
            // We also skip the next line;
            fileStream >> bypass >> bypass >> bypass >> bypass;

            // Since we are focusing on a single block, we just ignore the first numNodes lines
            // Because they are just numbers from 1 to numNodes
            for (int vI = 0; vI < numNodes; ++vI)
            {
                int vIndex;
                fileStream >> vIndex;
            }

            // Read vertex coordinates
            for (int vI = 0; vI < numNodes; ++vI)
            {
                fileStream >> TV(vI, 0) >> TV(vI, 1) >> TV(vI, 2);
            }
        }

        // Check for section indicating tetrahedra connectivity
        if (line.find("$Elements") != std::string::npos)
        {
            int bypass, numElements;
            fileStream >> bypass >> numElements >> bypass >> bypass; // Total nodes
            TT.resize(numElements, 4);
            // We also skip the next line;
            fileStream >> bypass >> bypass >> bypass >> bypass;

            // Read tetrahedra connectivity
            for (int elemI = 0; elemI < numElements; ++elemI)
            {
                int elemIndex;
                fileStream >> elemIndex >> TT(elemI, 0) >> TT(elemI, 1) >>
                    TT(elemI, 2) >> TT(elemI, 3);
            }
            // Convert indices to zero-based indexing
            TT.array() -= 1;
        }

        // Check for section indicating surface triangles
        if (line.find("$Surface") != std::string::npos)
        {
            int elemAmt;
            fileStream >> elemAmt; // Read the number of surface triangles
            SF.resize(elemAmt, 3);

            // Read surface triangles
            for (int triI = 0; triI < elemAmt; ++triI)
            {
                fileStream >> SF(triI, 0) >> SF(triI, 1) >> SF(triI, 2);
            }
            // Convert indices to zero-based indexing
            SF.array() -= 1;
        }
    }

    // If surface triangles are not provided in the file and flag is set to true, find surface triangles
    if (SF.rows() == 0)
    {
        std::cout << "Finding the surface triangle mesh for " << filePath << std::endl;
        igl::boundary_facets(TT, SF);
        // Convert indices to zero-based indexing
        SF.array() -= 1;
    }

    // Print the loaded mesh information
    std::cout << "tet mesh loaded with " << TV.rows() << " nodes, " << TT.rows() << " tets, and " << SF.rows() << " Surface triangles" << std::endl;

    return true;
}