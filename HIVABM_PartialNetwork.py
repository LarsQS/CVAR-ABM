#!/usr/bin/env python2.3
"""
Module that carries all definitions 
for the ABM of drug use and HIV. \n

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
from copy import deepcopy
import unittest
import numpy as np

try: from HIVABM_Population import PopulationClass
except ImportError:
    raise ImportError("Can't import PopulationClass")

try: import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("You must install matplotlib (http://matplotlib.sourceforge.net/)")
try: import networkx as nx
except ImportError:
    raise ImportError("You must install NetworkX (http://networkx.lanl.gov/) ")

class ScaleFreeNetworkClass(PopulationClass):

    def __init__(self, N = 10000, m_0 = 3):
        """
        :Purpose:
            This is the base class used to generate the social network 
            for the other agents, i.e. . The class inherits from the PopulationClass.

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
        SpecialAgents = set(self.IDU_agents).union(set(self.MSM_agents))
        SpecialAgents = SpecialAgents.union(set(self.NIDU_agents))
        NormalAgents = set(range(self.PopulationSize)).difference(SpecialAgents)
        NormalAgents = list(NormalAgents)
        self.NetworkSize = len(NormalAgents)
        # scale free Albert Barabsai Graph
        self.G = nx.barabasi_albert_graph(self.NetworkSize,m_0)


    def plot_DegreeDistribution(self, graph):
        """ 
        Plot the node degree distribution of the graph. \n 
        INPUT: networkX graph
        """
        print("\tNetwork size = "+str(self.NetworkSize))
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

    def get_Graph(self):
        """
        Return random assortative graph produced by ``set_assortative_graph``.
        """		
        return self.G

if __name__=='__main__':
    """unittest"""
    myNetworkObj = ScaleFreeNetworkClass(N=10000, m_0 = 3)
    G = myNetworkObj.get_Graph()
    myNetworkObj.plot_DegreeDistribution(G)