import java.util.LinkedList;
import java.util.Queue;

public class Path{
	
	//should be implemented with a LinkedList
	Stack<Integer> pathVertices;
	int pathWeight;

	public Path(Stack<Integer> pathVertices,int pathWeight){
		this.pathVertices = pathVertices;
		this.pathWeight = pathWeight;
	}	
	
	public void addPath(Stack<Integer> path, int weight){
		pathVertices.push(path);
		this.pathWeight += weight;
	}

	public void addVertex(int v, int weight){
		pathVertices.push(v);
		this.pathWeight += weight;
	}

	public void removeVertex(int weight){
		pathVertices.pop()
		this.pathWeight -= weight;
	}

	public int getWeight(){
		return this.pathWeight;
	}

	public Stack<Integer> getPath(){
		return this.pathVertices;
	}

}
