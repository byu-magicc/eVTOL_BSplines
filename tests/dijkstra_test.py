import numpy as np
from copy import deepcopy




#defines the function to get the minimum distance index from the uinitialized nodes
def getMinDistanceIndex(completeTable: list,
                        uninitializedNodes: list):

    minDistance = np.inf
    minDistanceNode = -1

    for tempList in completeTable:
        #gets the current distance and current node
        currentNode = tempList[0]
        currentDistance = tempList[1]
        previousNode = tempList[2]

        #if the current distance is less than the best found min distance, AND if
        #this node is in the uninitialized list
        if (currentDistance < minDistance) and (currentNode in uninitializedNodes):
            minDistance = currentDistance
            minDistanceNode = currentNode
    
    return minDistance, minDistanceNode



def main():
    #creates the list of possible combinations of nodes and their costs
    nodes = [0,1,2,3,4]

    node_connection_costs = {}
    node_connection_costs[(0,1)] = 5
    node_connection_costs[(0,2)] = 7
    node_connection_costs[(0,3)] = 12
    node_connection_costs[(1,4)] = 5
    node_connection_costs[(1,3)] = 9
    node_connection_costs[(2,3)] = 5
    node_connection_costs[(2,4)] = 4
    node_connection_costs[(3,4)] = 5


    #initializes the nodes section list
    nodes_sectionedList = []
    for i in range(len(nodes)):
        nodes_sectionedList.append([])

    #goes through and creates a list of the nodes each connects to and the cost
    for key, value in node_connection_costs.items():

        #gets the key first item
        key_start = key[0]
        key_end = key[1]

        #creates the new item
        newItem = (key_end, value)

        #gets the list within the sectioned list
        nodes_sectionedList[key_start].append(newItem)

        potato = 0


    #creates the lists of initialized and uninitialized nodes
    uninitializedNodes = deepcopy(nodes)
    initializedNodes = []

    #creates the main table of the connections back to zero
    mainTable = []
    for i in range(len(nodes)):
        #appends the index, the Shortest distance, and the previous Node
        #initializes the distance to zero, and the previous node to -1 if it's the first one
        if i == 0:
            tempTable = [i, 0, -1]
        #if it's the other ones, we initialize the shortest distance and previous nodes to infinity
        else:
            tempTable = [i, np.inf, np.inf]
        mainTable.append(tempTable)


    #iterates while len length of the length of the initialized nodes is not the length of the nodes
    while len(initializedNodes) != len(nodes):

        #gets the min distance node and its distance from the uninitialized nodes
        minDistance_cost, minDistanceNode = getMinDistanceIndex(completeTable=mainTable,
                                                           uninitializedNodes=uninitializedNodes)

        #now that we have the min distance node, we go through and update the nodes in the list from the 
        #things we have now.
        nodeConnections = nodes_sectionedList[minDistanceNode]
        #gets the current 

        for connection_pair in nodeConnections:

            #gets the end index
            endIndex = connection_pair[0]
            incremental_cost = connection_pair[1]

            #gets the total candidate cost
            candidate_cost = minDistance_cost + incremental_cost

            #gets the section of the main table
            mainTable_endIndex_section = mainTable[endIndex]

            #gets the current end index cost and parent
            current_endIndex_cost = mainTable_endIndex_section[1]
            
            #if our new candidate cost is less than the current end index cost,
            #we rewire
            if candidate_cost < current_endIndex_cost:
                #creates the new temp table
                newEntry = [endIndex, candidate_cost, minDistanceNode]
                #sets the main table with this entry
                mainTable[endIndex] = newEntry

            potato = 0

        #now that we've iterated through the connections from this node, we remove the min distance node from the 
        #uninitialized list and add it to the initialized list
        uninitializedNodes.remove(minDistanceNode)
        initializedNodes.append(minDistanceNode)

        tomato = 0


    niceCrispyBacon = 0

    #gets the end node index
    endNodeIndex = nodes[-1]
    startNodeIndex = nodes[0]

    currentIndex = deepcopy(endNodeIndex)
    #gets the total cost
    totalCost = (mainTable[endNodeIndex])[1]
    #creates the list of indices
    indices = [endNodeIndex]

    while currentIndex != startNodeIndex:
        
        #gets the current table section
        currentTableSection = mainTable[currentIndex]

        #gets the parent of the current index
        currentNode_parent = currentTableSection[2]
        #appends this to the list
        indices.append(currentNode_parent)
        #sets the current index as the parent
        currentIndex = currentNode_parent

    indices.reverse()


    po_ta_toes = 0



if __name__ == "__main__":
    main()