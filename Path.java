import java.util.Stack;

public class Path implements Cloneable{
	
	//should be implemented with a LinkedList
	Stack<Integer> pathVertices;
	double pathWeight;

	public Path(Stack<Integer> pathVertices,double pathWeight){
		this.pathVertices = pathVertices;
		this.pathWeight = pathWeight;
	}	
	
	public void addPath(Stack<Integer> path, double weight){
		pathVertices.addAll(path);
		this.pathWeight += weight;
	}

	public void addVertex(int v, double weight){
		pathVertices.push(v);
		this.pathWeight += weight;
	}

	public void removeVertex(double weight){
		pathVertices.pop();
		this.pathWeight -= weight;
	}

	public double getWeight(){
		return this.pathWeight;
	}

	public Stack<Integer> getPath(){
		return this.pathVertices;
	}

	public int pop(){
		return this.pathVertices.pop();
	}

	public int peek(){
		return this.pathVertices.peek();
	}

	public Path clone(){
		try{
			return (Path)super.clone();
		}
		catch(CloneNotSupportedException e)
		{
			throw new RuntimeException ("This class does not implement Cloneable");
		}
	}

}
