//Nicholas Guidice

#include <iostream>
#include <vector>
#include <queue>

using namespace std;


const int inf = 1e9; // A large value to represent infinity

//Structure to represent an edge
struct Edge {
    int vertex; //The destination vertex of the edge
    int weight; //The weight (cost) of the edge

    Edge(int v, int w) : vertex(v), weight(w) {}
};

//Class representing the graph using an adjacency list
class Graph {
public: 
    int numVertices; //Total number of vertices in the graph
    vector<vector<Edge>> adjList; //Adjacency list to store edges


    Graph(int vertices) : numVertices(vertices) {
        adjList.resize(vertices); //Resize the adjacency list to the number of vertices 
    }

    //Function to add an edge to the graph
    void addEdge(int u, int v, int weight) {
        adjList[u].emplace_back(v, weight); //Add an edge from vertex u to vertex v with the specified weight
    }
};

//Class for a priority queue using a binary heap
class PriorityQueue {
private:
    vector<pair<int, int>> heap; //Min-heap to store pairs of (distance, vertex)
    vector<int> verIndex;      //Index to keep track of vertices in the heap

public:
    PriorityQueue(int size) {
        verIndex.resize(size, -1); //Initialize all vertex indices to -1 (not in the heap)
    }

    //Function to push a vertex into the priority queue
    void push(int vertex, int distance) {
         //If the vertex is already in the queue, update its distance
        if (verIndex[vertex] != -1) {
            decreaseKey(vertex, distance);
            return;
        }

        //Otherwise add the new vertex to the heap
        heap.emplace_back(distance, vertex);
        verIndex[vertex] = heap.size() - 1; //Store the index of the vertex
        percolateUp(heap.size() - 1); //Maintain heap properties 
    } 

    //Function to pop the vertex with the minimum distance
    pair<int, int> pop() {
        if (heap.empty()) return make_pair(-1, -1); //Return a pair of -1 if the heap is empty
 
        pair<int, int> minElement = heap[0]; //Get the minimum element
        verIndex[minElement.second] = -1; //Mark the vertex as not in the heap
        heap[0] = heap.back(); //Move the last element to the root
        heap.pop_back();//Remove the last element

         // If the heap is not empty maintain the heap properties   
        if (!heap.empty()) {
            verIndex[heap[0].second] = 0;
            percolateDown(0);
        }
        return minElement; //Return the minimum element
    }

    //Check if the priority queue is empty
    bool isEmpty() const {
        return heap.empty();
    }

private:
//Percolate down the element at index index
    void percolateDown(int index) {
        int size = heap.size();
        while (true) {
            int left = 2 * index + 1; //Left child index
            int right = 2 * index + 2; //Right child index
            int smallest = index;
            //Check if the left child is smaller than the current smallest
            if (left < size && heap[left].first < heap[smallest].first)
                smallest = left;
            //Check if the right child is smaller than the current smallest
            if (right < size && heap[right].first < heap[smallest].first)
                smallest = right;
             //If the smallest is still the current index, the heap property is satisfied
            if (smallest == index) break;

            swap(heap[index], heap[smallest]); //Swap the current element with the smallest child
            verIndex[heap[index].second] = index;
            verIndex[heap[smallest].second] = smallest;
            index = smallest;
        }
    }

    //Percolate up the element at index index
    void percolateUp(int index) {
        while (index > 0) {
            int parent = (index - 1) / 2; //Parent Index
            if (heap[index].first >= heap[parent].first) break; //If current is greater than or equal to parent break
            swap(heap[index], heap[parent]); //Swap with parent
            verIndex[heap[index].second] = index;
            verIndex[heap[parent].second] = parent;
            index = parent; //Move up the heap
        }
    }

    //Decrease the key of a vertex
    void decreaseKey(int vertex, int newDistance) {
        int index = verIndex[vertex];
        if (index != -1 && newDistance < heap[index].first) {//If the vertex is in the heap and the new distance is smaller update it
            heap[index].first = newDistance;
            percolateUp(index);
        }
    }
};

//Dijkstra's algorithm implementation
pair<int, vector<int>> dijkstra(Graph& graph, int start, int end) {
    //Initialize distances with max values and previous vertices with -1
    vector<int> distances(graph.numVertices, inf);
    vector<int> previous(graph.numVertices, -1);
    distances[start] = 0; //Distance to the start vertex is 0

    PriorityQueue pq(graph.numVertices);
    pq.push(start, 0); //Push the start vertex into the queue

    while (!pq.isEmpty()) {
        
        pair<int, int> minElement = pq.pop(); 
        int currentDistance = minElement.first; 
        int currentVertex = minElement.second; 

        if (currentVertex == end) break; //Break if the end is reached 

         //Iterate over each neighbor of the current vertex
        for (size_t i = 0; i < graph.adjList[currentVertex].size(); ++i) {
            Edge edge = graph.adjList[currentVertex][i];  //Get the edge
            int neighbor = edge.vertex; //Neighbor vertex
            int weight = edge.weight; //Weight of the edge
            int distance = currentDistance + weight; //Calculate new distance

            // If the new distance is smaller update distances and previous vertex
            if (distance < distances[neighbor]) {
                distances[neighbor] = distance;
                previous[neighbor] = currentVertex;
                pq.push(neighbor, distance); //Push neighbor into the queue
            }
        }
    }


    return make_pair(distances[end], previous); // Return the shortest distance to the end vertex and the previous vertices
        
}

// Function to get the path from start to end
vector<int> getPath(int end, const vector<int>& previous) {
    vector<int> path; //Vector to hold the path
    for (int i = end; i != -1; i = previous[i]) {
        path.push_back(i); //Add the vertex to the path
    }
    vector<int> finalPath; 
    for (int i = path.size() - 1; i >= 0; --i) {
        finalPath.push_back(path[i]); // Add in reverse order
    }
    return finalPath;
}

//Main function
int main(int argc, char* argv[]) {
    //Check command line arguments are correct
    if (argc != 3) {
        cerr << "Invalid Arguments \n";
        return 1;
    }

    int startVertex = stoi(argv[1]);
    int endVertex = stoi(argv[2]);

    //Read the number of vertices and edges from input
    int numVertices, numEdges;
    cin >> numVertices >> numEdges;

    Graph graph(numVertices); //Create Graph with given amount of vertices
    
    //Read each edge and add it to the graph
    for (int i = 0; i < numEdges; ++i) {
        int u, v, weight;
        cin >> u >> v >> weight;
        graph.addEdge(u, v, weight);
    }

    //Call Dijkstra's
    pair<int, vector<int>> result = dijkstra(graph, startVertex, endVertex); 
    int distance = result.first; //Get the shortest distance
    vector<int> previous = result.second; //Get the previous vertices

    if (distance == inf) { //Check if the end vertex is reachable
        cout << "not connected\n";
    } else {
        cout << "Distance: " << distance << "\n"; 
        vector<int> path = getPath(endVertex, previous);
        cout << "Path: ";
        for (int vertex : path) {
            cout << vertex << " "; // Output the vertices in the path
        
    }
     cout << "\n";
    }

    return 0;
}
