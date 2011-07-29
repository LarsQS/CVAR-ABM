# NETWORKS
#
# Create network
#
#
# written by:
# Lars Seemann
# lseemann@uh.edu

import matplotlib.pyplot as plt
import matplotlib

from numpy import array, zeros
import networkx as nx


# Construct network model
class NetworkModel():
	def __init__(self, N = 1000):
		"""Initialize new network graph"""
		# test input
		if type(N) is not int:
			raise ValueError("Number of agents N must be integer")
		# lists
		self.PopulationSize = N
		self.graph = []
		self.adjacencyList = []
		
	def ErdosRenyi(self, p = 0.02):
		"""create Erdos-Renyi graph"""
		N = self.PopulationSize
		graph = nx.erdos_renyi_graph(n = N, p = p)
		self.graph = graph
		self.adjacencyList = graph.adjacency_list()
		
	def BarabasiAlbert(self, m = 4):
		"""create Barabasi Albert graph"""
		# test input
		if type(m) is not int:
			raise ValueError("Input m must be integer")	
		# main procedure
		N = self.PopulationSize
		graph = nx.barabasi_albert_graph(n = N, m = m)
		self.graph = graph
		self.adjacencyList = graph.adjacency_list()
		
	def WattsStrogatz(self, nei = 4, p = 0.2):
		"""create Watts Strogatz graph"""
		# test input
		if type(nei) is not int:
			raise ValueError("Number of neigbors nei must be integer")
		# main procedure
		N = self.PopulationSize
		graph = nx.watts_strogatz_graph(n = N, k = nei, p = p)
		self.graph = graph
		self.adjacencyList = graph.adjacency_list()
		
	def Lattice(self, k = 5):
		"""
		CREATE LATTICE
		
		each agent is part of a kxk neighborhood
		within each neighborhood agents are organized on a lattice and have exactly 4 nearest neighbors
		Periodic (circular) boundary conditions, i.e. agents on the edge are connected 'to the other site'
		This behavior can be changes in Graph.Lattice function: circular=False
		
		"""
		# Test input
		if type(k) is not int:
			raise ValueError("Size k of neighborhood lattice must be integer")
		N = self.PopulationSize
		NumNeiH = N/(k*k)
		if N%(k*k) !=0:
			raise ValueError("Number of agents must be dividable into neighborhood block")
		else:
			print ("Number of neighborhoods: " + str(NumNeiH))
		# MAIN PROCEDURE
		# Idea: create one kxk lattice, use it again with different agent indices
		# Creat single lattice
		g = nx.grid_2d_graph(k, k, periodic=True)
		# convert node label (kxk) into integer (just count through)
		gLattice = nx.convert.convert_node_labels_to_integers(g)
		adjListLattice = gLattice.adjacency_list()
		# Create NumNeiH neighborhoods within agents are organized in lattice
		# Idea: use above lattice for all neighborhoodblock, just change indices in the lattice
		CompleteAdjList = []									# complete Adj List of all agents
		for i in range(NumNeiH):								# for each neighborhood
			for iAgent in range(k*k):							# count through all agents in one neighborhood
				tmpListneiboors = []							# temporary list of neighbors of that one agent
				for agent in adjListLattice[iAgent]:			# get all neighbors of that one agent
					tmpListneiboors.append(agent + (i*(k*k)) )	# get neighbor list bu add i*(k*k) to each iAgnet neighbor
				CompleteAdjList.append(tmpListneiboors)
		# get Link List from adjacency List
		LinkList = []
		for iAgent in range(len(CompleteAdjList)):
			for agent in CompleteAdjList[iAgent]:
				t = iAgent, agent
				LinkList.append(t)
		# graph from list of links
		graph = nx.Graph()
		graph.add_edges_from(LinkList)
		self.graph = graph
		self.adjacencyList = graph.adjacency_list()
		#print self.adjacencyList
					
	def PlotGraph(self):
		"""plot the graph"""
		plt.figure(figsize=(8,8))
        # with nodes sized by degree
		# node_color=[float(H.degree(v)) for v in H]
		# layout:
		#pos=nx.spring_layout(self.graph)
		pos = nx.circular_layout(self.graph)
		#pos=nx.shell_layout(self.graph)
		#pos=nx.spectral_layout(self.graph)
		
		node_color = 'b'
		edge_color = 'k'
		node_shape = 'o'
		NodeSize=[]
		for v in self.graph:
			NodeSize.append((10*self.graph.degree(int(v)))**(1.0))
		
		# draw:
		nx.draw(self.graph, pos, 
			node_size = NodeSize, 
			node_color = node_color,
			node_shape = node_shape,
			edge_color = edge_color,
			alpha=0.5, with_labels=False)
		#nx.draw_networkx_nodes(self.graph,pos,node_size=NodeSize)
		#nx.draw_networkx_edges(self.graph,pos,alpha=0.4)
		plt.axis('equal')
		plt.axis('off')
		plt.show()
		
	def PlotDegreeDistribution(self):
		""" plot node degree distribution """
		print(" Populationsize = "+str(self.PopulationSize))
		degreeList = array( self.graph.degree().values() )	
		degreemax = (max(degreeList)+1)
		degreehist = zeros(degreemax)
		x_degree = range(degreemax)

		for n in x_degree:
			degreehist[n] = len(degreeList[degreeList == n])

		# delete zeros
		ix = degreehist != 0
		degreehist = degreehist[ix]
		x_degree = array(x_degree)
		x_degree = x_degree[ix]
			
		# print node degree distribution	
		plt.loglog(x_degree, degreehist, 'bo')
		plt.title('Node degree distribution')
		plt.ylabel('Frequency')
		plt.xlabel('Node degree')
		plt.axis('equal')
		plt.axis([min(degreeList)-(0.05*max(degreeList)), max(degreeList)+(0.05*max(degreeList)), 0.9, max(degreehist)+(0.05* max(degreehist))])
		plt.show()
		
	def get_adjlist(self):
		""" return adjacency list"""
		return self.adjacencyList

			
if __name__=='__main__':
	
	N = 10
	myModelNetwork = NetworkModel(N = N)
	#myModelNetwork.Lattice(k=5)
	#myModelNetwork.BarabasiAlbert(m=2)
	myModelNetwork.WattsStrogatz(nei = 4, p = 0.2)
	
	print 'Test output:'
	#print myModelNetwork.PopulationSize
	#print myModelNetwork.adjacencyList
	print myModelNetwork.graph
	
	myModelNetwork.PlotGraph()
	myModelNetwork.PlotDegreeDistribution()
	