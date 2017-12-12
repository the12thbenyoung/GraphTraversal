
//Modified from Graph.java from the package EDU.colorado.graphs

import java.util.LinkedList;
import java.util.Queue;

public class Graph implements Cloneable {
	/* Adjacency matrix:
	 * - for edges[i][j]: i is source vertex, j is target vertex
	 * - each element in edges is the path weight between the current i and j
	 * - if edges[i][j] = -1, then no path exists between i and j
	 */
	private double[][] edges;
	
	//Labels array: each source vertex i has a label contained in labels[i]
	private Object[] labels;

	public Graph(double[][] adjacencies){
		edges = adjacencies;
		labels = new Object[adjacencies.length]; //all values initially null
	}
	
	// **********EDGES**********
	
	public boolean isEdge(int source, int target){
		return (edges[source][target] != -1);
	}

	// **********WEIGHTS**********
	
	public double getWeight(int source, int target){
		return edges[source][target];
	}
	
	// **********LABELS**********
	
	public void setLabel(int vertex, Object newLabel){
		labels[vertex] = newLabel;
	} 
	
	public Object getLabel(int vertex){
		return labels[vertex];
	}
	
	public int size(){
		return labels.length;
	}
	
	// **********CLONING**********
	
	public Object clone() {
	   Graph answer;
	   try {
		   answer = (Graph)super.clone();
	   }
	   catch(CloneNotSupportedException e) {  
		   // This would indicate an internal error in the Java runtime system
		   // because super.clone always works for an Object.
	      throw new InternalError(e.toString());
	   }
	   
	   answer.edges = (double[][]) edges.clone();
	   answer.labels = (Object[]) labels.clone();
	   
	   return answer;
	}
	
	// **********DEPTH FIRST SEARCHING**********
	
	public static void depthFirstPrint(Graph g, int start) {
		// Boolean array for whether or not vertex has been marked (traversed once)
		boolean[] marked = new boolean [g.size()];
		depthFirstRecurse(g, start, marked);  
	}
	
	public static void depthFirstRecurse(Graph g, int v, boolean[] marked){
		
		// Starting vertex:
		marked[v] = true; // Mark
		System.out.println(g.getLabel(v)); // Process
		
		int[] connections = g.neighbors(v);
		
		// Traverse all reighboring vertices
		for (int i = 0; i < connections.length; i++){
			int nextNeighbor = connections[i];
			// Check if reighbor vertex is marked
			if (!marked[nextNeighbor]){
				depthFirstRecurse(g, nextNeighbor, marked);
			}
		} 
	}
	
	// **********BREADTH FIRST SEARCHING**********
	
	public static void breadthFirstPrint(Graph g, int start){
		// Boolean array for whether or not vertex has been marked (traversed once)
		boolean[] marked = new boolean[g.size()];
		breadthFirst(g, start, marked);
	}
	
	public static void breadthFirst(Graph g, int v, boolean[] marked){
		// Queue to store unmarked neighbor vertices to process
		Queue<Integer> q = new LinkedList<Integer>();
		
		// Starting vertex:
		System.out.println(g.getLabel(v)); // Process
		marked[v] = true; // Mark
		q.add(v); // Place in queue
		
		// Run until queue is empty (no more unmarked vertices in scope)
		while(!q.isEmpty()){
			// Remove vertex from front
			// Traverse all neighboring vertices
			int[] connections = g.neighbors(q.remove());
			for (int i = 0; i < connections.length; i++) {
				int nextNeighbor = connections[i];
			    if (!marked[nextNeighbor]) { // Check if neighbor vertex is marked
					System.out.println(g.getLabel(nextNeighbor));// Process
			    	marked[nextNeighbor] = true; // Mark
			    	q.add(nextNeighbor); // Place in queue
				}
			}
		} 
	}
	
	// **********NEIGHBORING VERTICES**********
	
	public int[] neighbors(int vertex){
	   int i;
	   int count;
	   int[] answer;
	   
	   // First count how many edges have the vertex as their source
	   count = 0;
	   for (i = 0; i < labels.length; i++){
	      if (edges[vertex][i] != -1)
	         count++;
	   }
	        
	   // Allocate the array for the answer
	   answer = new int[count];
	   
	   // Fill the array for the answer
	   count = 0;
	   for (i = 0; i < labels.length; i++){
	      if (edges[vertex][i] != -1)
	         answer[count++] = i;
	   }
	   return answer;
	}
}