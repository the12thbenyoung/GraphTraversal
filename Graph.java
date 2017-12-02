import java.util.LinkedList;
import java.util.Queue;

//File: Graph.java from the package EDU.colorado.graphs
//Complete documentation is available from the Graph link in
//http://www.cs.colorado.edu/~main/docs/

//package edu.colorado.graphs;

public class Graph implements Cloneable
{
// Invariant of the Graph class:
//   1. The vertex numbers range from 0 to labels.length-1.
//   2. For each vertex number i, labels[i] contains the label for vertex i.
//   3. For any two vertices i and j, edges[i][j] is true if there is a
//      vertex from i to j; otherwise edges[i][j] is false.  
private double[][] edges;
private Object[ ] labels;

public Graph(int n){
   edges = new double[n][n];  // All values initially -1
   for (int i = 0; i < n; i++) {
	   for (int j = 0; j < n; j++) {
		   edges[i][j] = -1;
	   }
   }
   labels = new Object[n];     // All values initially null
}

public void addEdge(int source, int target, double weight){
   edges[source][target] = weight;
}

public double returnWeight(int source, int target){
	return edges[source][target];
}


public Object clone( )
{  // Clone a Graph object.
   Graph answer;
   
   try
   {
      answer = (Graph) super.clone( );
   }
   catch (CloneNotSupportedException e)
   {  // This would indicate an internal error in the Java runtime system
      // because super.clone always works for an Object.
      throw new InternalError(e.toString( ));
   }
   
   answer.edges = (double [ ][ ]) edges.clone( );
   answer.labels = (Object [ ]) labels.clone( );
   
   return answer;
}


public static void depthFirstPrint(Graph g, int start)
{
   boolean[ ] marked = new boolean [g.size( )];
   
   depthFirstRecurse(g, start, marked);
   
}


public static void depthFirstRecurse(Graph g, int v, boolean[ ] marked)
{
   int[ ] connections = g.neighbors(v);
   int i;
   int nextNeighbor;
   
   marked[v] = true;
   System.out.println(g.getLabel(v));
   
   // Traverse all the neighbors, looking for unmarked vertices:
   for (i = 0; i < connections.length; i++)
   {
      nextNeighbor = connections[i];
      if (!marked[nextNeighbor])
         depthFirstRecurse(g, nextNeighbor, marked);
   } 
}

public static void breadthFirstPrint(Graph g, int start) {
	boolean[] marked = new boolean[g.size()];
	breadthFirst(g, start, marked);
}

public static void breadthFirst(Graph g, int v, boolean[] marked) {
	Queue<Integer> q = new LinkedList<Integer>();
	
	//starting vertex
	System.out.println(g.getLabel(v)); //process
	marked[v] = true; //mark
	q.add(v); //place in queue
	
	//remaining vertices
	while(!q.isEmpty()){
		//remove vertex from front
		int[] connections = g.neighbors(q.remove()); //find neighbors of vertex
		for (int i = 0; i < connections.length; i++) {
			int nextNeighbor = connections[i];
		    if (!marked[nextNeighbor]) { //for each unmarked neighbor
				System.out.println(g.getLabel(nextNeighbor));//process
		    	marked[nextNeighbor] = true; //mark
		    	q.add(nextNeighbor); //place in queue
			}
		}
	} 
}

public Object getLabel(int vertex)
{
   return labels[vertex];
}


public boolean isEdge(int source, int target)
{
   return (edges[source][target] != -1);
}

public int[ ] neighbors(int vertex)
{
   int i;
   int count;
   int[ ] answer;
   
   // First count how many edges have the vertex as their source
   count = 0;
   for (i = 0; i < labels.length; i++)
   {
      if (edges[vertex][i] != -1)
         count++;
   }
        
   // Allocate the array for the answer
   answer = new int[count];
   
   // Fill the array for the answer
   count = 0;
   for (i = 0; i < labels.length; i++)
   {
      if (edges[vertex][i] != -1)
         answer[count++] = i;
   }
   
   return answer;
}
           

public void removeEdge(int source, int target)   
{
   edges[source][target] = -1;
}

public void setLabel(int vertex, Object newLabel)
{
   labels[vertex] = newLabel;
}   

public int size( )
{
   return labels.length;
}
     
}

