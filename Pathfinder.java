import java.util.Arrays;
import java.util.Stack;
public class Pathfinder{
	public static void main(String args[]){
		Graph g;	
		Path p;
		long[] runTimes;
		int numTrials;
		int numVerticesStart;
	        int numVerticesEnd;
		int numVerticesStep;	

		//vary edges for random graph-------------------------------------------------------
/*		numTrials = 10000;
 *		int avgEdgesStart = 3;
*	        int avgEdgesEnd = 8;
*		int avgEdgesStep = 1;	
*		int numVertices = 10;
*		System.out.println("Brute force algorithm for random graph with numVertices = " + numVertices);		
*		
*		long startTime,endTime;
*		runTimes = new long[(avgEdgesEnd - avgEdgesStart)/avgEdgesStep + 1];
*
*		System.out.print("Average Edges");
*		for(int avgEdges = avgEdgesStart; avgEdges <= avgEdgesEnd; avgEdges += avgEdgesStep){
*			System.out.print("\t" + avgEdges);
*			g = randomGraph(numVertices,avgEdges,1000);
*			startTime = System.currentTimeMillis();
*			for(int i = 0; i < numTrials; i++){
*				if(i%(numTrials/10) == 0){  //makes 10 new graphs through the trial to average out runtimes
*					g = randomGraph(numVertices,avgEdges,1000);
*				}
*				p = bruteForceShortestPath(g);
*			}
*			endTime = System.currentTimeMillis();
			runTimes[(avgEdges - avdEdgesStart)/avgEdgesStep] = endTime - startTime;
*		}
*		System.out.print("\nRun Time(ms)");
*		for(int i = 0; i < runTimes.length; i++){
*			System.out.print("\t" + runTimes[i]);	
*		}
*		System.out.println();
*/		
		//vary number of vertices for random graph-------------------------------------------------------
/*		numTrials = 10000;
 *		numVerticesStart = 7;
*	        numVerticesEnd = 12;
*		numVerticesStep = 1;	
*		int avgEdges = 6;
*		System.out.println("Brute force algorithm for random graph with avgEdges = " + avgEdges);
*		
*		long startTime,endTime;
*		runTimes = new long[(numVerticesEnd - numVerticesStart)/numVerticesStep + 1];
*
*		System.out.print("Number of Vertices");
*		for(int numVertices = numVerticesStart; numVertices <= numVerticesEnd; numVertices += numVerticesStep){
*			System.out.print("\t" + numVertices);
*			g = randomGraph(numVertices,avgEdges,1000);
*			startTime = System.currentTimeMillis();
*			for(int i = 0; i < numTrials; i++){
*				if(i%(numTrials/10) == 0){  //makes 10 new graphs through the trial to average out runtimes
*					g = randomGraph(numVertices,avgEdges,1000);
*				}
*				p = bruteForceShortestPath(g);
*			}
*			endTime = System.currentTimeMillis();
			runTimes[(numVertices - numVerticesStart)/numVerticesStep] = endTime - startTime;
*		}
*		System.out.print("\nRun Time(ms)");
*		for(int i = 0; i < runTimes.length; i++){
*			System.out.print("\t" + runTimes[i]);	
*		}
*		System.out.println();
*/	
	
		//vary number of vertices for linear graph-------------------------------------------------------
/*		numVerticesStart = 100;
*	        numVerticesEnd = 700;
*		numVerticesStep = 100;	
*		numTrials = 10000;
*		System.out.println("Brute force algorithm for linear graph");
*		
*		long startTime,endTime;
*		runTimes = new long[(numVerticesEnd - numVerticesStart)/numVerticesStep + 1];
*
*		System.out.print("Number of Vertices");
*		for(int numVertices = numVerticesStart; numVertices <= numVerticesEnd; numVertices += numVerticesStep){
*			System.out.print("\t" + numVertices);
*			g = linearGraph(numVertices);
*			startTime = System.currentTimeMillis();
*			for(int i = 0; i < numTrials; i++)
*				p = bruteForceShortestPath(g);
*			endTime = System.currentTimeMillis();
*			runTimes[(numVertices - numVerticesStart)/numVerticesStep] = endTime - startTime;
*		}
*		System.out.print("\nRun Time(ms)");
*		for(int i = 0; i < runTimes.length; i++){
*			System.out.print("\t" + runTimes[i]);	
*		}
*
*		System.out.println();
*/
		//vary dead-end edges (neighbors) per vertex for branched graph-------------------------------------------------------
/*		numTrials = 1;
*		int neighborsStart = 1;
*		int neighborsEnd = 10;
*		int neighborsStep = 1;	
*		int numMainVertices = 100;
*		System.out.println("Brute force algorithm for branching graph with numMainVertices = " + numMainVertices);		
*		
*		long startTime,endTime;
*		runTimes = new long[(neighborsEnd - neighborsStart)/neighborsStep + 1];
*
*		System.out.print("Dead-end Edges");
*		for(int neighbors = neighborsStart; neighbors <= neighborsEnd; neighbors += neighborsStep){
*			System.out.print("\t" + neighbors);
*			g = branchingGraph(numMainVertices,neighbors);
*			startTime = System.currentTimeMillis();
*			for(int i = 0; i < numTrials; i++){
*				p = bruteForceShortestPath(g);
*			}
*			endTime = System.currentTimeMillis();
*			runTimes[(neighbors - neighborsStart)/neighborsStep] = endTime - startTime;
*		}
*		System.out.print("\nRun Time(ms)");
*		for(int i = 0; i < runTimes.length; i++){
*			System.out.print("\t" + runTimes[i]);	
*		}
*		System.out.println();
*
*/			
		//vary number of main-path vertices for branching graph-------------------------------------------------------
		numTrials = 10000;
		numVerticesStart = 100;
                numVerticesEnd = 700;
		numVerticesStep = 100;	
		int neighbors = 3;
		System.out.println("Brute force algorithm for branching  graph with neighbors per main vertex = " + neighbors);
		
		long startTime,endTime;
		runTimes = new long[(numVerticesEnd - numVerticesStart)/numVerticesStep + 1];

		System.out.print("Number of Vertices");
		for(int numVertices = numVerticesStart; numVertices <= numVerticesEnd; numVertices += numVerticesStep){
			System.out.print("\t" + numVertices);
			g = branchingGraph(numVertices,neighbors);
			startTime = System.currentTimeMillis();
			for(int i = 0; i < numTrials; i++){
				p = bruteForceShortestPath(g);
			}
			endTime = System.currentTimeMillis();
			runTimes[(numVertices - numVerticesStart)/numVerticesStep] = endTime - startTime;
		}
		System.out.print("\nRun Time(ms)");
		for(int i = 0; i < runTimes.length; i++){
			System.out.print("\t" + runTimes[i]);	
		}
		System.out.println();
		
		}
	

	/*linearGraph
	 *Creates a linear graph with only one path from start to end 
	 *@Parameters:
	 *size - the number of vertices in the path
	 */
	public static Graph linearGraph(int size) {
		double[][] matrix = new double[size][size];
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (j == i+1) matrix[i][j] = 1;
				else matrix[i][j] = -1;
			}
		}
		return new Graph(matrix);
	}
	
	/*branchingGraph
	 *Creates a branching graph - similar to a linear graph except each vertex in the main path has 
	 *lower weight edges leading to dead ends grafted onto it. Designed to confuse Dijkstra's algorithm
	 *and decrease its efficiency
	 *@Parameters:
	 *size - the number of vertices in the main chain (the path from start to end vertex)
	 *neighbors - the number of vertices grafted onto each vertex in the main chain
	 */
	public static Graph branchingGraph(int size, int neighbors) {
		double[][] matrix = new double[size + neighbors*(size-1)+1][size + neighbors*(size-1)+1];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix.length; j++) {
				matrix[i][j] = -1;
			}
		}
		//first <size> vertices in graph are the main path, size + (i-1)*neighbors + 1 in the index of the main path vertex
		//i's first neighbor
		for (int i = 0; i < size; i++) {
			for (int j = i; j < size; j++) {
				if (j == i+1) {
					matrix[i][j] = 1;  //sets edges in main path to larger weight
					for(int n = 0; n < neighbors; n++){
						matrix[i][size + (i)*neighbors + 1 + n] = 0;  //sets edges to dead-end vertices to smaller weight
					}                                                                                                     
				}
			}
		}
		return new Graph(matrix);
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
		//fills in last row with -1s (end node shouldn't have any paths leaving from it)
		for(int col = 0; col < size; col++){
			edges[size-1][col] = -1;
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
		
		return g;
	}


public static int minDistance(Graph g, double[] distance, boolean[] marked)
    {
        // Initialize min value
        double min = Integer.MAX_VALUE;
        int minIndex = -1;
 
        for (int i = 0; i < g.size(); i++)
            if (distance[i] != -1 && marked[i] == false && distance[i] <= min) {
                min = distance[i];
                minIndex = i;
            }
 
        return minIndex;
    }
 
    // A utility function to print the constructed distance array
    public static void printSolution(Graph g, double[] distance, int end)
    {
        System.out.println("Vertex\tDistance from Source");
        for (int i = 0; i <= end; i++)
            System.out.println(i+"\t"+distance[i]);
    }
 
    // Function that implements Dijkstra's single source shortest path
    // algorithm for a graph represented using adjacency matrix
    // representation
    public static void dijkstra(Graph g, int start, int end)
    {
    	double[][] matrix = new double[g.size()][g.size()];
    	for (int i = 0; i < matrix.length; i++) {
    		for (int j = 0; j < matrix.length; j++) {
    			matrix[i][j] = g.getWeight(i, j);
    		}
    	}
    	
        double distance[] = new double[g.size()]; // The output array. distance[i] will hold
                                 // the shortest distance from start to i
 
        // marked[i] will true if vertex i is included in shortest
        // path tree or shortest distance from start to i is finalized
        boolean[] marked = new boolean[g.size()];
 
        // Initialize all distances as INFINITE and stpSet[] as false
        for (int i = 0; i < g.size(); i++)
        {
            distance[i] = Integer.MAX_VALUE;
            marked[i] = false;
        }
 
        // distance of source vertex from itself is always 0
        distance[start] = 0;
 
        // Find shortest path for only the vertices
        for (int i = 0; i < end; i++)
        {
        		// Pick the minimum distance vertex from the set of vertices
	            // not yet processed. u is always equal to start in first
	            // iteration.
	            int u = minDistance(g, distance, marked);
	 
	            // Mark the picked vertex as processed
	            marked[u] = true;
	 
	            // Update distance value of the adjacent vertices of the
	            // picked vertex.
	            for (int v = 0; v < g.size(); v++) {
	                // Update distance[v] only if is not in marked, there is an
	                // edge from u to v, and total weight of path from start to
	                // v through u is smaller than current value of distance[v]
	                if (!marked[v] && matrix[u][v]!=-1 && distance[u] != Integer.MAX_VALUE && distance[u]+matrix[u][v] < distance[v])
	                    	distance[v] = distance[u] + matrix[u][v];
	            }
        	}
 
        // print the constructed distance array
        printSolution(g, distance, end);
    }





	public static Path bruteForceShortestPath(Graph g){
		//marked makes sure that this path doesn't visit any vertex it already visited
		boolean[] marked = new boolean[g.size()];

		return bruteForceRecurse(g, 0, marked);  
	}
	
	public static Path bruteForceRecurse(Graph g, int v, boolean[] marked){
		Stack<Integer> s = new Stack<Integer>();
		s.push(v);
		if(v == marked.length - 1){ //if this vertex is the target (target is always last node)
			return new Path(s,0);
		}


		marked[v] = true;
		double maxWeight = (double)Integer.MAX_VALUE;
		int[] connections = g.neighbors(v);

		Path shortestPath = new Path(s,maxWeight);
		Path tempPath;
		boolean isDeadEnd = true;

		// Traverse all neighboring vertices
		for (int i = 0; i < connections.length; i++){
			int nextNeighbor = connections[i];
			// Check if neighbor vertex is marked
			if (!marked[nextNeighbor]){
				isDeadEnd = false;
				tempPath = bruteForceRecurse(g, nextNeighbor, marked);
				//if this new path is shorter than the current shortest one, make it the new shortest one
				if(tempPath.getWeight() + g.getWeight(v,nextNeighbor) < maxWeight){
					maxWeight = tempPath.getWeight() + g.getWeight(v,nextNeighbor);
					shortestPath = tempPath;
				}
			}

		} 

		//after all edges leaving from this vertex have been searched, remove it from marked
		marked[v] = false;
		
		//add this vertex to the shortest path
		if(isDeadEnd) //if this is a dead end...
			shortestPath.addVertex(v,(double)Integer.MAX_VALUE); //the path weight should be effectively infinite
		else{
			int nextVertex = shortestPath.peek();
			shortestPath.addVertex(v,g.getWeight(v,nextVertex)); //adds the current vertex and the weight of the edge between it and the first vertex in the shortest path to the end
		}
		return shortestPath; 
	}
}	
