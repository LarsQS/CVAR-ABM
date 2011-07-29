#!/usr/bin/env python2.3
"""
Module that carries the defintions for the Population for 
the ABM of drug use and HIV. \n

Author:		Lars Seemann \n
Email:		lseemann@uh.edu \n 
Date:		2011-01-28 \n

Copyright (c) 2010, under the Simplified BSD License. \n
For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php \n
All rights reserved.
"""
__author__="Lars Seemann (lseemann@uh.edu)"

import random
import copy
import numpy as np
from copy import deepcopy
import unittest

try: from HIVABM_Population import PopulationClass
except ImportError:
    raise ImportError("Can't import PopulationClass")

try: import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("You must install matplotlib (http://matplotlib.sourceforge.net/)")
try: import networkx as nx
except ImportError:
    raise ImportError("You must install NetworkX (http://networkx.lanl.gov/) ")


class  SocialNetworkClass(PopulationClass):
    """
    :Purpose:
    	This class inherits from the class :py:class:`PopulationClass`. It therefore 
    	inherits the population and embeds all agents in a random social network. 
        
    :Input:	
    	N : int
    		Number of agents. Default: 10000
                
    	m_0: int	
    		Number of nodes each node is connected to in preferential
    		attachment step
                
    	:py:class:`PopulationClass` : Inherited	
        
    :Attributes:
    	:py:attr:`G` : networkX graph object
    	:py:attr:`pool_of_similar_nodes` \n
    	:py:attr:`repeated_nodes` \n
    	:py:attr:`repeated_nodes_IDU` \n
    	:py:attr:`repeated_nodes_NIDU` \n
    	:py:attr:`repeated_nodes_ND` \n
    	:py:attr:`repeated_nodes_HIV` \n
    	All attributes from :py:class:`PopulationClass`
        
    :Methods:
    	All methods from :py:class:`PopulationClass` \n
    	:py:meth:`_initialize_network` \n
    	:py:meth:`_set_assortative_graph` \n
    	:py:meth:`get_assortative_graph` \n
    	:py:meth:`visualize_network` \n
    	:py:meth:`plot_DegreeDistribution` \n
    	:py:meth:`get_AdjacencyList` \n
    """
    def __init__(self, N = 10000, m_0 = 3):
        """
        :Purpose:
          This is the base class used to generate the social network 
          and neighborhoods. The class inherits from the PopulationClass.

        :Input:	
            N : int
              Number of agents. Default: 10000
            m_0: int	
              Number of nodes each node is connected to in preferential
              attachment step
        """
        if type(N) is not int:			
            raise ValueError(('Population size must be integer,\
					n = %s, not %s')%(string(N), type(N)))	
        else: pass
        if m_0 not in range(10):
            raise ValueError('m_0 must be integer smaller than 10')
        else: self.m_0 = m_0

        PopulationClass.__init__(self, n = N)	# Create population
        self.repeated_nodes = {'All':[], 'IDU':[], 'NIDU':[], 'ND':[], 'HIV':[]}
        self.G = nx.Graph()
        self._set_assortative_graph()		# create graph

    def _initialize_network_old(self):
        """
        :Purpose:
        	Initialize random network generator with a network consisting of:
        	- two IDU agent
        	- two NIDU agent
        	- two ND agent
        	and connect the three agents completely (4links each).
        :Output:
        	G : networkX object
        	Seed network
        """
        self.repeated_nodes = {'All':[], 'IDU':[], 'NIDU':[], 'ND':[], 'HIV':[]}

        targets = []
        targets.append(random.choice(self.ND_agents))
        targets.append(random.choice(self.IDU_agents))
        targets.append(random.choice(self.NIDU_agents))
        targets.append(random.choice(self.MSM_agents))
        targets.append(random.choice(self.HIV_agents))
        targets.append(random.choice(self.AIDS_agents))

        # The network seed is going to be fully connected
        # therefore we have as many nodes as the *targets*
        # is long (n), and we have n-1 edges coming from
        # each node ( total nr of edges = n(n-1)/2 )
        m = len(targets) - 1

        for agent in targets:
            self.repeated_nodes['All'].extend([agent]*m)
            DrugType = self.get_agent_charactersitic(agent,'Drug Type')
            self.repeated_nodes[DrugType].extend([agent]*m)
            if agent in self.HIV_agents: self.repeated_nodes['HIV'].extend([agent]*m)

        G=nx.Graph()			# Completely connect seed network 
        G.name="assortative_graph(%s)"%(self.PopulationSize)
        G.add_nodes_from(targets)
        NumTargets = len(targets)
        for u in xrange(NumTargets):
            for v in xrange(u+1,NumTargets):
                G.add_edge(targets[u],targets[v])
        return G

    def _initialize_network(self, m, pool):
        """
        :Purpose:
          Initialize random network generator with a network consisting of:
          - two IDU agent
          - two NIDU agent
          - two ND agent
          and connect the three agents completely (4links each).
        :Output:
          G : networkX object
        Seed network
        """
        targets = []
        while len(targets)<m:
            targets.append(random.choice(pool))
        # The network seed is going to be fully connected
        # therefore we have as many nodes as the *targets*
        # is long (n), and we have n-1 edges coming from
        # each node ( total nr of edges = n(n-1)/2 )
        mm = len(targets) - 1
        for agent in targets:
            self.repeated_nodes['All'].extend([agent]*mm)
        G=nx.Graph()	# Completely connect seed network 
        G.name="assortative_graph(%s)"%(self.PopulationSize)
        G.add_nodes_from(targets)
        NumTargets = len(targets)
        for u in xrange(NumTargets):
            for v in xrange(u+1,NumTargets):
                G.add_edge(targets[u],targets[v])
        self.G = deepcopy(G)
        if len(self.G)!=m:
            raise ValueError("Check initialization! len(self.G) = %d, m = %d"%
                             (len(self.G), m))

    def _random_subset(self,seq,m):
        """
        :Purpose:
        Return m unique elements from seq.\n
        This differs from random.sample which can return repeated
        elements if seq holds repeated elements.
        :Input:	
        seq : list
        List from which the m random elements are chosen.
        m : int	
        Number of elemts.
        :Output:
        targets : set
        """
        return set(random.sample(set(seq), m))

    def _set_assortative_graph(self):
        """
        Updated 05/01/2011
        """
        p_IDU = 0.8					    # Probability of IDUs to attach to HIV
        p_MSM = 0.8					    # Probability of MSMs to attach to HIV
        m_0_std = 2                                     # Nr of added edges at each step
        m_0_msm = range(10)                             # Nr of added edges at each step
        m_0_idu = range(26,41)                          # Nr of added edges at each step
        self.repeated_nodes = {'All':[], 'IDU':[], 'MSM':[]}
        self.repeated_nodes['MSM'].extend(self.MSM_agents)
        self.repeated_nodes['IDU'].extend(self.IDU_agents)

        # Barabasi Albert for non MSM, IDU agents.MSM and IDU agents handled seperately.
        NormalAgents = list(set(range(self.PopulationSize)).difference(
            set(self.IDU_agents).union(set(self.MSM_agents))))
        self._initialize_network(m_0_std, NormalAgents)
        # Exclude agents that were used to initialize the network:
        for agent in self.G.nodes():
            if agent in NormalAgents: NormalAgents.remove(agent)
            else: raise ValueError("agent %s not in NormalAgents"%agent)
        # Add 'normal' agents to network
        for source in NormalAgents:
            targets = self._random_subset(self.repeated_nodes['All'], m_0_std)
            self.G.add_edges_from(zip([source]*m_0_std,targets))    # Add edges from source to m nodes
            self.repeated_nodes['All'].extend([source]*m_0_std)	# Add source node to lists
            self.repeated_nodes['All'].extend(targets)		# Add target nodes to lists
            for ag in targets:
                if self.get_agent_charactersitic(ag,'Sex Type') == 'MSM':
                    self.repeated_nodes['MSM'].append(ag)
                if self.get_agent_charactersitic(ag,'Drug Type') == 'IDU':
                    self.repeated_nodes['IDU'].append(ag)
        if len(self.G)!=(len(NormalAgents)+2):
            raise ValueError("len(G)=%d,Num Normal Agents=%d"%(len(self.G), (len(NormalAgents)+2)))

        # MSM agents excluding IDU
        MSMAgents = list(set(self.MSM_agents).difference(set(self.IDU_agents)))
        for source in MSMAgents:
            m = random.choice(m_0_msm)
            if random.random() < p_MSM:
                targets = self._random_subset(self.repeated_nodes['MSM'], m)
            else: targets = self._random_subset(self.repeated_nodes['All'], m)
            # Bookkeeping
            self.G.add_edges_from(zip([source]*m_0_std,targets))    # Add edges from source to m nodes
            self.repeated_nodes['All'].extend([source]*m_0_std)	# Add source node to list
            self.repeated_nodes['MSM'].extend([source]*m_0_std)	# Add source node to MSM list
            self.repeated_nodes['All'].extend(targets)		# Add target nodes to lists
            for ag in targets: 
                if self.get_agent_charactersitic(ag,'Sex Type') == 'MSM':
                    self.repeated_nodes['MSM'].append(ag)
                if self.get_agent_charactersitic(ag,'Drug Type') == 'IDU':
                    self.repeated_nodes['IDU'].append(ag)
            if len(self.G)!=(len(NormalAgents)+2+len(MSMAgents)):
                raise ValueError("len(G)=%d, Expected:%d"%(
                    len(self.G), (len(NormalAgents)+2+len(MSMAgents))))
        # IDU agents
        for source in self.IDU_agents:
            m = random.choice(m_0_idu)
            if random.random() < p_IDU:
                targets = self._random_subset(self.repeated_nodes['IDU'], m)
            else: targets = self._random_subset(self.repeated_nodes['All'], m)
            # Bookkeeping
            self.G.add_edges_from(zip([source]*m_0_std,targets))    # Add edges from source to m nodes
            self.repeated_nodes['All'].extend([source]*m_0_std)	# Add source node to list
            self.repeated_nodes['IDU'].extend([source]*m_0_std)	# Add source node to IDU list
            self.repeated_nodes['All'].extend(targets)		# Add target nodes to lists
            for ag in targets: 
                if self.get_agent_charactersitic(ag,'Sex Type') == 'MSM':
                    self.repeated_nodes['MSM'].append(ag)
                if self.get_agent_charactersitic(ag,'Drug Type') == 'IDU':
                    self.repeated_nodes['IDU'].append(ag)
            if len(self.G)!=(len(NormalAgents)+2+len(MSMAgents)+len(self.IDU_agents)):
                raise ValueError("len(G)=%d, Expected:%d"%(len(self.G),
                                                           (len(NormalAgents)+2+len(MSMAgents)+len(self.IDU_agents))))
        # Check consistency:
        if len(self.G)!=self.PopulationSize:
            NormalAgents = list(set(range(self.PopulationSize)).difference(
                set(self.IDU_agents).union(set(self.MSM_agents))))
            print "Network size: %d"%len(self.G)
            print "PopulationSize: %d"%self.PopulationSize
            print "Normal Agents: %d"%len(NormalAgents)
            print "MSM Agents: %d"%len(MSMAgents)
            print "IDU Agents: %d"%len(self.IDU_agents)
            raise ValueError("len(G)=%d, Population size = %d"%(
                len(self.G), self.PopulationSize))


    def _set_assortative_graph_old(self, p=0.95):
        """
        :Purpose:
        	Create random graph using a modified Barabasi-Albert preferential 
        	attachment model. In the control model agents are placed in the network 
        	according to demographic assignment with preferential grouping among 
        	agents with identical drug usage.
        :Input:
        	N :	int	
        		Number of agents
        	m_0: int	
        		Number of nodes each node is connected to in preferential
        		attachment step
        :Details:
        	A graph of n nodes is grown by attaching new nodes each with m
        	edges that are preferentially attached to existing nodes with high
        	degree. 
        	The assortative behavior is modeled by randomly attaching m nodes to the new
        	incoming node. This attachment step is biased so that nodes with the same 
        	characteristics have a higher chance to be attached to the new graph. 
        	This biased preferential attachment step causes assortative mixing 
        	which reflects social segregation. 
        :References:
        	.. [1] A. L. Barabasi and R. Albert "Emergence of scaling in
        	   random networks", Science 286, pp 509-512, 1999.
        """
        PreferentialProbability = p
        p_IDU = 0.8					# Probability of IDUs to attach to HIV
        p_MSM = 0.4					# Probability of MSMs to attach to HIV
        m_0_std = 1                                     # Nr of added edges at each step
        m_0_msm = range(10)                             # Nr of added edges at each step
        m_0_idu = range(26-41)                          # Nr of added edges at each step
        allAgents = range(self.PopulationSize)		# List of all agents
        self.repeated_nodes = {'All':[], 'IDU':[], 'NIDU':[], 'ND':[], 'HIV':[]}
        self.G = self._initialize_network()		# Seed network


        # Remove initial agents from pool of possible connections
        # because they are already connected
        for agent in self.G.nodes():
            if agent in allAgents: allAgents.remove(agent)
            else: raise ValueError("agent %s not in allAgents"%agent)
        for source in allAgents:			# Start adding the other n-m_0 nodes to network G
            DrugType = self.get_agent_charactersitic(source,'Drug Type')
            SexType = self.get_agent_charactersitic(source,'Sex Type')
            if DrugType == 'IDU': m_0 = m_0_idu
            elif SexType == 'MSM': m_0 = m_0_msm
            else: m_0 = m_0_std
            targets = set()			        # Choose new targets
            while len(targets)<m_0:
                tmp_rnd = random.random()	        # Only vertex degree considered
                if len(self.G.nodes()) < 10 or tmp_rnd < (1-PreferentialProbability):
                    x=random.choice(self.repeated_nodes['All'])
                elif DrugType == 'IDU':
                    if tmp_rnd < (1-PreferentialProbability) + p_IDU: 
                        x=random.choice(self.repeated_nodes['HIV'])
                    else: x=random.choice(self.repeated_nodes['All'])
                elif SexType == 'MSM':
                    if tmp_rnd < (1-PreferentialProbability) + p_MSM: 
                        x=random.choice(self.repeated_nodes['HIV'])
                    else: x=random.choice(self.repeated_nodes['All'])
                else: x=random.choice(self.repeated_nodes['All'])
                if type(x) != int: raise ValueError("Target must be int! x=%s"%x)
                targets.add(x)
            self.G.add_edges_from(zip([source]*m_0_std,targets))    # Add edges from source to m nodes
            ## Update Lists
            self.repeated_nodes['All'].extend([source]*m_0_std)	    # Add source node to lists
            self.repeated_nodes['All'].extend(targets)		    # Add target nodes to lists
            self.repeated_nodes[DrugType].extend([source]*m_0_std)
            if source in self.HIV_agents: self.repeated_nodes['HIV'].extend([source]*m_0_std) 
            for ag in targets:
                DrugType_ag = self.get_agent_charactersitic(ag,'Drug Type')
                self.repeated_nodes[DrugType_ag].append(ag)
                if ag in self.HIV_agents: self.repeated_nodes['HIV'].append(ag)			

    def get_assortative_graph(self):
        """
        Return random assortative graph produced by ``set_assortative_graph``.
        """		
        return self.G

    def visualize_network(self, graph):
        """
        Visualize the network using the spring layout (default). \n
        INPUT: networkX graph
        """
        print("Plotting...")
        plt.figure(figsize=(8,8))
        # with nodes sized by degree
        # node_color=[float(H.degree(v)) for v in H]
        # layout:
        pos=nx.spring_layout(graph)
        #pos = nx.circular_layout(graph)
        #pos=nx.shell_layout(graph)
        #pos=nx.spectral_layout(graph)

        edge_color = 'k'
        node_shape = 'o'
        # node color to indicate ethnicity, shape occupation
        node_color = []
        for v in graph:
            if v in self.IDU_agents:
                node_color.append('r')
            elif v in self.NIDU_agents:
                node_color.append('b')
            elif v in self.ND_agents:
                node_color.append('g')
            else:
                raise ValueError("Agent ID invalid!")

        # node size indicating popularity
        NodeSize=[]
        for v in graph:
            NodeSize.append((10*graph.degree(int(v)))**(1.0))

        # draw:
        nx.draw(graph, pos, 
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

    def plot_DegreeDistribution(self, graph):
        """ 
        Plot the node degree distribution of the graph. \n 
        INPUT: networkX graph
        """
        print(" Populationsize = "+str(self.PopulationSize))
        degreeList = np.array( graph.degree().values() )	
        x_degree = np.arange(max(degreeList)+1)	# include 0 and max
        degreehist = np.bincount(degreeList)
        ix = degreehist != 0					# delete zeros
        degreehist = degreehist[ix]
        x_degree = x_degree[ix]

        # print node degree distribution	
        plt.loglog(x_degree, degreehist, 'bo')
        plt.title('Node degree distribution')
        plt.ylabel('Frequency')
        plt.xlabel('Node degree')
        plt.axis('equal')
        plt.axis([min(degreeList)-(0.05*max(degreeList)), 
                  max(degreeList)+(0.05*max(degreeList)),
                  0.9, max(degreehist)+(0.05* max(degreehist))])
        plt.show()

    def get_AdjacencyList(self):
        """ Return the adjacency list of the graph. """
        return self.G.adjacency_list()

    def plot_stratified_DegreeDistribution(self):
        """ 
        :Purpose:
        	Plot the node degree distribution of the graph. 
        """
        graph = self.G
        degreeList_IDU = np.array( graph.degree(self.IDU_agents).values() )
        degreeList_MSM = np.array( graph.degree(self.MSM_agents).values() )
        Other_agents = list(set(self.Agents.keys()) - (set(self.IDU_agents) | set(self.MSM_agents)) )
        degreeList_Other = np.array( graph.degree(Other_agents).values() )

        x_degree_IDU = np.arange(max(degreeList_IDU)+1)	# include 0 and max
        degreehist_IDU = np.bincount(degreeList_IDU)
        degreehist_IDU = degreehist_IDU/float(np.sum(degreeList_IDU))
        ix = degreehist_IDU != 0			# delete zeros
        degreehist_IDU = degreehist_IDU[ix]
        x_degree_IDU = x_degree_IDU[ix]

        x_degree_MSM = np.arange(max(degreeList_MSM)+1)	# include 0 and max
        degreehist_MSM = np.bincount(degreeList_MSM)
        degreehist_MSM = degreehist_MSM/float(np.sum(degreeList_MSM))
        ix = degreehist_MSM != 0			# delete zeros
        degreehist_MSM = degreehist_MSM[ix]
        x_degree_MSM = x_degree_MSM[ix]

        x_degree_Other = np.arange(max(degreeList_Other)+1)	# include 0 and max
        degreehist_Other = np.bincount(degreeList_Other)
        degreehist_Other = degreehist_Other/float(np.sum(degreeList_Other))
        ix = degreehist_Other != 0			        # delete zeros
        degreehist_Other = degreehist_Other[ix]
        x_degree_Other = x_degree_Other[ix]


        # print node degree distribution	
        #plt.loglog(x_degree_IDU, degreehist_IDU, '--go', label='IDU')
        #plt.loglog(x_degree_MSM, degreehist_MSM, '--bs', label='MSM')
        #plt.loglog(x_degree_Other, degreehist_Other, '--r^', label='Other')
        plt.plot(x_degree_IDU, degreehist_IDU, '--go', label='IDU')
        plt.plot(x_degree_MSM, degreehist_MSM, '--bs', label='MSM')
        plt.plot(x_degree_Other, degreehist_Other, '--r^', label='Other')
        plt.legend()
        plt.xlim((0,40))
        plt.title('Stratified Node Degree Distribution')
        plt.ylabel('Frequency')
        plt.xlabel('Node degree / Number of contacts')
        plt.show()

    def assess_risk_contacts(self):
        """ 
        :Purpose:
        	Assess the probability for risk contacts for 
        	MSM, IDU and all other agents. Plot the results.
        """
        Count = {'IDU':[] , 'MSM':[] , 'Other':[]}
        for agent in self.Agents.keys():
            SexType = self.get_agent_charactersitic(agent, 'Sex Type')
            DrugType = self.get_agent_charactersitic(agent, 'Drug Type')
            if SexType == 'MSM':
                tmpCount = 0
                for partner in self.G[agent]:
                    if self.get_agent_charactersitic(partner, 'HIV') != 0:
                        tmpCount += 1
                Count['MSM'].append(tmpCount)
            if DrugType == 'IDU':
                tmpCount = 0
                for partner in self.G[agent]:
                    if self.get_agent_charactersitic(partner, 'HIV') != 0:
                        tmpCount += 1
                Count['IDU'].append(tmpCount)
            if SexType != 'MSM' and DrugType != 'IDU':
                tmpCount = 0
                for partner in self.G[agent]:
                    if self.get_agent_charactersitic(partner, 'HIV') != 0:
                        tmpCount += 1
                Count['Other'].append(tmpCount)	
        return Count
        plt.plot(np.bincount(Count['IDU']), '--go', label='IDU')
        plt.plot(np.bincount(Count['MSM']), '--bs', label='MSM')
        plt.plot(np.bincount(Count['Other']), '--r^', label='Other')
        plt.legend()


class TestClassMethods(unittest.TestCase):
    """ 
    :Purpose:
    	Unittest for testing the methods of the HIVABM network class
    """
    def setUp(self):
        """
        :Purpose:
        	Tests that all models from setup pass inspection. 
        	``setUp`` is perfomed before each method.
        """
        self.N_pop = 20000

    def test_NetworkInitialization(self):
        """ Test: Testing the population"""
        print "\t __Testing the network initialization"	
        myNetwork = SocialNetworkClass(N=self.N_pop)
        G = myNetwork._initialize_network()
        NumIDU = NumNIDU = NumND = 0
        for node in G.nodes():
            if node in myNetwork.IDU_agents: NumIDU +=1
            elif node in myNetwork.NIDU_agents: NumNIDU +=1
            elif node in myNetwork.ND_agents: NumND +=1
            else: raise ValueError("Node undefined!")
        self.assertTrue(G.number_of_edges() == 15)

def assess_sampling():
    """ 
    :Purpose:
    	Test the sampling for MSM, IDU and all other agents
    """
    N_mc = 10
    TimeMax = 10
    Npop = 10000
    # Store results in list for each Monte Carlo run
    Count_allcontacts = {'IDU': [] , 'MSM': [] , 'Other': []}
    Count_HIVcontacts = {'IDU': [] , 'MSM': [] , 'Other': []}
    for i in range(N_mc):
        if (i+1)%(N_mc/10)==0: print "\t%d %%"%int(N_mc/(i+1))
        # Count absolute number of contacts over lifetime
        tmpCount_allcontacts = {'IDU': 0 , 'MSM': 0 , 'Other': 0}
        tmpCount_HIVcontacts = {'IDU': 0 , 'MSM': 0 , 'Other': 0}
        myNetwork = SocialNetworkClass(N=Npop)
        for t in range(TimeMax):
            for agent in myNetwork.Agents.keys():
                SexType_agent = myNetwork.get_agent_charactersitic(agent, 'Sex Type')
                DrugType_agent = myNetwork.get_agent_charactersitic(agent, 'Drug Type')
                HIVType_agent = myNetwork.get_agent_charactersitic(agent, 'HIV')
                partner = random.choice(myNetwork.G[agent].keys())
                SexType_partner = myNetwork.get_agent_charactersitic(partner, 'Sex Type')
                DrugType_partner = myNetwork.get_agent_charactersitic(partner, 'Drug Type')
                HIVType_partner = myNetwork.get_agent_charactersitic(partner, 'HIV')
                # Consider agent -> partner
                if SexType_agent == 'MSM':
                    tmpCount_allcontacts['MSM'] += 1
                    if HIVType_partner != 0: tmpCount_HIVcontacts['MSM'] += 1
                if DrugType_agent == 'IDU':
                    tmpCount_allcontacts['IDU'] += 1
                    if HIVType_partner != 0: tmpCount_HIVcontacts['IDU'] += 1
                if SexType_agent != 'MSM' and DrugType_agent != 'IDU':
                    tmpCount_allcontacts['Other'] += 1
                    if HIVType_partner != 0: tmpCount_HIVcontacts['Other'] += 1
                # Consider partner -> agent
                if SexType_partner == 'MSM':
                    tmpCount_allcontacts['MSM'] += 1
                    if HIVType_agent != 0: tmpCount_HIVcontacts['MSM'] += 1
                if DrugType_partner == 'IDU':
                    tmpCount_allcontacts['IDU'] += 1
                    if HIVType_agent != 0: tmpCount_HIVcontacts['IDU'] += 1
                if SexType_partner != 'MSM' and DrugType_partner != 'IDU':
                    tmpCount_allcontacts['Other'] += 1
                    if HIVType_agent != 0: tmpCount_HIVcontacts['Other'] += 1
        for type in ['IDU', 'MSM', 'Other']:
            Count_allcontacts[type].append(tmpCount_allcontacts[type])
            Count_HIVcontacts[type].append(tmpCount_HIVcontacts[type])

    plt.plot(np.bincount(Count_HIVcontacts['IDU']), '--go', label='IDU')
    plt.plot(np.bincount(Count_HIVcontacts['MSM']), '--bs', label='MSM')
    plt.plot(np.bincount(Count_HIVcontacts['Other']), '--r^', label='Other')
    plt.legend()


if __name__=='__main__':
    """unittest"""
    #unittest.main()
    assess_sampling()
