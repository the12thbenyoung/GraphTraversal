import java.util.ArrayList;
import java.util.Arrays;
import java.util.PriorityQueue;

public class DijkstraGraph {
	Double[][] matrix;
	static int size;
	
	public DijkstraGraph(Double[][] m) {
		this.matrix = m;
		DijkstraGraph.size = m.length;
	}
	
	public static DijkstraGraph LinearGraph(int size) {
		Double[][] matrix = new Double[size][size];
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (j == i+1) matrix[i][j] = 1.0;
				else matrix[i][j] = null;
			}
		}
		//printMatrix(matrix);
		return new DijkstraGraph(matrix);
	}
	
	public static DijkstraGraph BranchingGraph(int size, int neighbors) {
		Double[][] matrix = new Double[size+neighbors*(size-1)][size+neighbors*(size-1)];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix.length; j++) {
				matrix[i][j] = null;
			}
		}
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (j == i+1) {
					matrix[i][j] = 1.0;
					for (int n = 0; n < neighbors; n++) {
						matrix[i][size + i*neighbors + n] = 0.0;
					}
				}
			}
		}
		
		//printMatrix(matrix);
		return new DijkstraGraph(matrix);
	}
	
	public static DijkstraGraph RandomGraph(int size, int averageNeighborSize, int maxEdgeWeight) {
		Double[][] matrix = new Double[size][size];
		for (int row = 0; row < size-1; row ++) {
			for (int col = 0; col < size; col++) {
				//gives vertex <row> a edge to vertex <col> avgEdgesPerVertex/size percent of the time
				if ((Math.random() < (double) averageNeighborSize / size) && row != col) {
					//edge weight is a random positive int less than maxEdgeWeight
					matrix[row][col] = (double)(int)(Math.random()*maxEdgeWeight);
				}
				else matrix[row][col] = null;
			}
		}
		//fills in last row with nulls (end node shouldn't have any paths leaving from it)
		for (int col = 0; col < size; col++)
			matrix[size-1][col] = null;

		//checks how many edges lead to last (end of path) vertex
		int edgesToLastVertex = 0;
		for (int row = 0; row < size - 1; row++) {
			if (matrix[row][size-1] != null)
				edgesToLastVertex++;
		}

		//if last vertex doesn't have at least 2 edges leading to it, add extra edges to it
		while (edgesToLastVertex < 2) {
			//if the second to last vertex doesn't already have an edge to the last vertex, give it one
			if (matrix[size-2][size-1] == null)
				matrix[size-2][size-1] = (double)(int)(Math.random()*maxEdgeWeight);
			//otherwise go to the vertex before that
			else matrix[size-3][size-1] = (double)(int)(Math.random()*maxEdgeWeight);
			edgesToLastVertex++;
		}
		//checks how many edges come from first (start of path) vertex
		int edgesFromFirstVertex = 0;
		for (int col = 0; col < size - 1; col++) {
			if (matrix[0][col] != null)
				edgesFromFirstVertex++;
		}

		//if first vertex doesn't have at least 2 edges leaving it, add extra edges
		while(edgesFromFirstVertex < 2) {
			//if the second vertex doesn't already have an edge from the first vertex, give it one
			if(matrix[0][1] == null)
			        matrix[0][1] = (double)(int)(Math.random()*maxEdgeWeight);
			//otherwise go to the vertex after that
			else matrix[0][2] = (double)(int)(Math.random()*maxEdgeWeight);
			edgesFromFirstVertex++;
		}
		//printMatrix(matrix);
		return new DijkstraGraph(matrix);
	}
		
	public void dijkstra(int endVertex) {
		boolean[] marked = new boolean[size];
		runDijkstra(endVertex, marked);
	}
	
	public void runDijkstra(int endVertex, boolean[] marked) {
		Double[] path = new Double[size];	// initially null: infinity
		Integer[] previous = new Integer[size];
		path[0] = 0.0;
		marked[0] = true;
		PriorityQueue<PathVertex> pq = new PriorityQueue<PathVertex>();
		pq.offer(new PathVertex(0, 0.0));
		while (!pq.isEmpty()) {
			PathVertex current = pq.poll();
			int vertex = current.vertex;
			double pathWeight = current.pathWeight;
			if (pathWeight == path[vertex]) {
				for (Integer i: getNeighbors(vertex)) {
					if (!marked[i]) {
						relaxEdge(path, previous, vertex, i);
						pq.add(new PathVertex(i, path[i]));
						marked[i] = true;
					}
				}
			}
		}
	}
	
	public void relaxEdge(Double[] path, Integer[] previous,
						  int sourceVertex, int targetVertex) {
		double weight = getEdgeWeight(sourceVertex, targetVertex);
		Double value = path[targetVertex];
		/* if distance is currently at infinity */
		if (value == null) {
			path[targetVertex] = path[sourceVertex] + weight;
			previous[targetVertex] = sourceVertex;
		}
		else if (value > path[sourceVertex] + weight) {
			path[targetVertex] = path[sourceVertex] + weight;
			previous[targetVertex] = sourceVertex;
		}
	}
	
	public ArrayList<Integer> getNeighbors(int vertex) {
		ArrayList<Integer> neighbors = new ArrayList<Integer>();
		for (int i = 0; i < size; i++) {
			if (getEdgeWeight(vertex, i) != null) {
				neighbors.add(i);
			}
		}
		return neighbors;
	}
	
	public Double getEdgeWeight(int sourceVertex, int targetVertex) {
		return matrix[sourceVertex][targetVertex];
	}
		
	public static void printMatrix(Double[][] matrix) {
		System.out.println();
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix.length; j++) {
				System.out.print(matrix[i][j] + "|");
			}
			System.out.println();
		}
		System.out.println();
	}
	
	public static void main(String[] args) {
		
		//Double[][] matrix = {{ null,	0.0,	null,	2.0,	4.0},
			//				   {null,	null,	null, 	null,	0.0},
				//			   {7.0,	4.0,	null, 	3.0,	null},
					//		   {9.0,	4.0,	5.0, 	null,	null},
						//	   {null,	null,	null,	null,	null}};
		//DijkstraGraph graph = new DijkstraGraph(matrix);
		//graph.dijkstra(4);
		//DijkstraGraph linear1 = LinearGraph(5);
		//linear1.dijkstra(4);
		//DijkstraGraph branching1 = BranchingGraph(3, 2);
		//branching1.dijkstra(2);
		//DijkstraGraph random1 = RandomGraph(10000, 3, 10);
		//random1.dijkstra(9999);
	}
	
	static class PathVertex implements Comparable<DijkstraGraph.PathVertex> {
		Integer vertex;
		Double pathWeight;
		public PathVertex(int vertex, Double path) {
			this.vertex = vertex;
			this.pathWeight = path;
		}
		public int compareTo(PathVertex other) {
			return this.pathWeight.compareTo(other.pathWeight);
		}
	}
}