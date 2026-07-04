#include "mesh.hpp"

#include "constants.hpp"
#include "utils.hpp"

void Mesh::insertNodesAtPositions(const std::vector<Vec2D> nodePos)
{
    const unsigned int lastIndexOfExistingNodes = nodes.nodesPos.size() - 1;
    const unsigned int numAddedNodes = nodePos.size();
    nodes.numNodes += numAddedNodes;
    nodes.nodesPos.insert(nodes.nodesPos.begin(), nodePos.begin(), nodePos.end());
    nodes.constraints.numFree += numAddedNodes;
    std::vector<unsigned int> newFreeIndices;
    newFreeIndices.reserve(numAddedNodes);
    for (unsigned int i = 0; i < numAddedNodes; i++) newFreeIndices.push_back(i+1+lastIndexOfExistingNodes);

    nodes.constraints.freeIndices.insert(nodes.constraints.freeIndices.begin(), newFreeIndices.begin(), newFreeIndices.end());

}
