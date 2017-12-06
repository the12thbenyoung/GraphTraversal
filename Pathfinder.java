import java.util.Arrays;
import java.util.Stack;
public class Pathfinder{
	public static void main(String args[]){
		Graph g = randomGraph(5,4,100);
		Path p = bruteForceShortestPath(g);
		Stack<Integer> s = p.getPath();
		int stackSize = s.size();
		for(int i = 0; i < stackSize; i++){
			System.out.println(s.pop());
		}
	}

	public static Path bruteForceShortestPath(Graph g){
		//marked makes sure that this path doesn't visit any vertex it already visited
		boolean[] marked = new boolean[g.size()];

		//path keeps track of which vertices are in the current path. thisPath[g.size()] holds the pathLength
		Stack<Integer> s = new Stack<Integer>();
		s.push(0);
		Path path = new Path(s,0);

		return bruteForceRecurse(g, 0, path, marked);  
	}
	
	public static Path bruteForceRecurse(Graph g, int v, Path path, boolean[] marked){
		//this edge is the edge between last vertex in the path and v
		Path thisPath = path.clone();

		if(v == marked.length - 1){ //if this vertex is the target (target is always last node)
			thisPath.addVertex(v,0);
			return thisPath;
		}


		double maxWeight = (double)Integer.MAX_VALUE;
		Path shortestPath = new Path(null,maxWeight);
		marked[v] = true;
		int[] connections = g.neighbors(v);

		Path tempPath;
		boolean isDeadEnd = true;

		// Traverse all neighboring vertices
		for (int i = 0; i < connections.length; i++){
			int nextNeighbor = connections[i];
			// Check if neighbor vertex is marked
			if (!marked[nextNeighbor]){
				isDeadEnd = false;
				tempPath = bruteForceRecurse(g, nextNeighbor, thisPath, marked);
				if(tempPath.getWeight() < maxWeight){
					maxWeight = tempPath.getWeight();
					shortestPath = tempPath;
				}
			}

		} 

		//after all edges leaving from this vertex have been searched, remove it from marked
		marked[v] = false;
		
		//add this vertex to the shortest path
		if(isDeadEnd) //if this is a dead end...
			shortestPath.addVertex(v,(double)Integer.MAX_VALUE); //the path weight should be effectively infinite
		else
			shortestPath.addVertex(v,g.getWeight(v,shortestPath.peek()));
		return shortestPath; 
	}	
	


	/*randomGraph - generates a pseudo-random graph
	 *parameters:
	 *size - the number of vertices in the graph
	 *avgEdgesPerVertex - average number of edges leaving each vertex
	 *maxEdgeWeight - maximum weight of any edge in the graph
	 *returns: a Graph object
	 */
	public static Graph randomGraph(int size,int avgEdgesPerVertex, int maxEdgeWeight){
		double[][] edges = new double[size][size];
		for(int row = 0; row < size-1; row ++){
			for(int col = 0; col < size; col++){
				//gives vertex <row> a edge to vertex <col> avgEdgesPerVertex/size percent of the time
				if((Math.random() < (double)avgEdgesPerVertex/size) && row != col){
					//edge weight is a random positive int less than maxEdgeWeight
					edges[row][col] = (int)(Math.random()*maxEdgeWeight);
				}
				else{
					edges[row][col] = -1;	
				}
			}
		}

		//checks how many edges lead to last (end of path) vertex
		int edgesToLastVertex = 0;
		for(int row = 0; row < size - 1; row++){
			if(edges[row][size-1] > 0)
				edgesToLastVertex++;
		}

		//if last vertex doesn't have at least 2 edges leading to it, add extra edges to it
		while(edgesToLastVertex < 2){
			//if the second to last vertex doesn't already have an edge to the last vertex, give it one
			if(edges[size-2][size-1] == -1)
			        edges[size-2][size-1] = (int)(Math.random()*maxEdgeWeight);
			//otherwise go to the vertex before that
			else
			        edges[size-3][size-1] = (int)(Math.random()*maxEdgeWeight);

			edgesToLastVertex++;
		}
		//checks how many edges come from first (start of path) vertex
		int edgesFromFirstVertex = 0;
		for(int col = 0; col < size - 1; col++){
			if(edges[0][col] > 0)
				edgesFromFirstVertex++;
		}

		//if first vertex doesn't have at least 2 edges leaving it, add extra edges
		while(edgesFromFirstVertex < 2){
			//if the second vertex doesn't already have an edge from the first vertex, give it one
			if(edges[0][1] == -1)
			        edges[0][1] = (int)(Math.random()*maxEdgeWeight);
			//otherwise go to the vertex after that
			else
			        edges[0][2] = (int)(Math.random()*maxEdgeWeight);

			edgesFromFirstVertex++;
		}


		Graph g = new Graph(edges);

		for(int i = 0; i < size; i++){
			g.setLabel(i,"V"+i);
		}
		
		for(int i = 0; i < size; i++){
			for(int j = 0; j < size; j++){
				System.out.print(edges[i][j] + "\t");
			}
			System.out.println();
		}
		
		return g;
	}

}
