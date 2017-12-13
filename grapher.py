
import numpy
import matplotlib.pyplot as p

#Random graph changing edge density ------------------------------------------------
def randGraphEdge(fig,ax):
    bruteForceTimes = [431,840,3709,9114,26562,107436]
    dijkstraTimes = [151,108,69,54,61,63]
    x = [3,4,5,6,7,8]

    ax.scatter(x,bruteForceTimes,label = "brute force algorithm")
    ax.scatter(x,dijkstraTimes,label = "Dijkstra's algorithm")

    ax.set(xlabel='Average number of edges per vertex', ylabel='Run time (ms)',
            title='Random Graph with 10 Vertices')

#Random graph changing number of vertices ------------------------------------------------
def randGraphVertex(fig,ax):
    bruteForceTimes = [942,1935,3656,7796,24736,80942]

    dijkstraTimes = [124,86,83,71,60,86]
    x = [7,8,9,10,11,12]
    ax.scatter(x,bruteForceTimes,label = "brute force algorithm")
    ax.scatter(x,dijkstraTimes,label = "Dijkstra's algorithm")
    ax.set(xlabel='Number of Vertices', ylabel='Run time (ms)',
            title='Random Graph with Average of 6 Edges per Vertex')


#Linear graph changing number of vertices ------------------------------------------------
def linearGraph(fig,ax):
    bruteForceTimes = [768,2023,4357,7827,11090,16242,21177]
    dijkstraTimes = [1422,4357,9643,17294,27144,38807,52520]
    x = [100,200,300,400,500,600,700]
    ax.scatter(x,bruteForceTimes,label = "brute force algorithm")
    ax.scatter(x,dijkstraTimes,label = "Dijkstra's algorithm")

    ax.set(xlabel='Number of Vertices', ylabel='Run time (ms)',
            title='Linear Graph')

    #Branching graph changing edge density ------------------------------------------------
def branchedGraphEdge(fig,ax):
    bruteForceTimes = [2098,3878,6736,9783,14240,18840,24460,29936,36537,43646]
    dijkstraTimes = [4061,8750,14394,22639,32527,43870,57049,71561,88164,107026]
    edges = [1,2,3,4,5,6,7,8,9,10]

    ax.scatter(edges,bruteForceTimes,label = "brute force algorithm")
    ax.scatter(edges,dijkstraTimes,label = "Dijkstra's algorithm")

    ax.set(xlabel='Number of Dead-end edges per main vertex', ylabel='Run time (ms)',
            title='Branched Graph with 100 main vertices')

    #Random graph changing edge density ------------------------------------------------
def branchedGraphVertex(fig,ax):
    x = [100,200,300,400,500,600,700]
    dijkstraTimes = [1443,5784,12913,22761,35378,50718,68978]
    bruteForceTimes = [791,2460,5210,9100,14030,19805,27058]

    ax.scatter(x,bruteForceTimes,label = "brute force algorithm")
    ax.scatter(x,dijkstraTimes,label = "Dijkstra's algorithm")

    ax.set(xlabel='Number of main-path vertices', ylabel='Run time (ms)',
            title='Branched Graph with 3 dead-end edges per main vertex')

    
fig,ax= p.subplots()
randGraphEdge(fig,ax)
#randGraphVertex(fig,ax)
#linearGraph(fig,ax)
#branchedGraphEdge(fig,ax)
#branchedGraphVertex(fig,ax)
ax.legend()
p.show()
