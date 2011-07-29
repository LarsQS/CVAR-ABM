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

class PopulationClass():
    """
    :Purpose:
        This class constructs and represents the model population

    :Input:	
    
        n : int
            Number of agents. Default: 10000

    :Attributes:
    
        :py:attr:`PopulationSize` : int
            Size of the population.

        :py:attr:`propIDU` : float
            Prevalence of intravenous drug users.

        :py:attr:`numIDU` : int
            Number of intravenous drug users.

        :py:attr:`propNIDU` : float
            Prevalence of non-intravenous drug users.

        :py:attr:`numNIDU` : int
            Number of non-intravenous drug users.

        :py:attr:`propND` : float
            Prevalence of not drug users.

        :py:attr:`numND` : int
            Number of not drug users.

        :py:attr:`Agents`: dict
            Dictionary of agents and their characteristics. The agents are
            the `keys' and a dicotionary of 'characteristic:value'
            pair is the entry.

        :py:attr:`IDU_agents`: list
            IDU drug users

        :py:attr:`NIDU_agents`: list
            NIDU drug users

        :py:attr:`ND_agents`: list
            ND drug users

        :py:attr:`MSM_agents`: list
            MSM agents

        :py:attr:`HIV_agents`: list
            HIV+ users

        :py:attr:`AIDS_agents`: list
            Users with AIDS

    :Methods:
        :py:meth:`get_agent_charactersitic` \n
        :py:meth:`_set_population` \n
        :py:meth:`get_agents` \n
        :py:meth:`get_info_DrugUserType` \n
        :py:meth:`get_info_HIV_IDU` \n
        :py:meth:`get_info_DrugSexType`
    """

    def __init__(self, n = 10000):
        """
        :Purpose:
            Initialize PopulationClass object.
        """
        if type(n) is not int:
            raise ValueError("Population size must be integer")
        else:
            self.PopulationSize = n
        # Parameters
        # drug user prevalence (proportion)
        self.propIDU  = 191.0/10000.0
        self.numIDU = int(self.propIDU * self.PopulationSize)
        self.propNIDU = 647.0/10000.0
        self.numNIDU = int(self.propNIDU * self.PopulationSize)
        self.propND = 9162.0/10000.0
        self.numND = self.PopulationSize - self.numIDU - self.numNIDU
        self.Agents = dict()	# Main Dictionary Agents
        # List of IDU, NIDU, NDs
        # shuffle all Agents
        allAgents = range(self.PopulationSize)
        random.shuffle(allAgents)
        # print ('Number of agents '+str(len(allAgents)))
        self.IDU_agents = deepcopy(allAgents[0:self.numIDU])
        self.NIDU_agents = deepcopy(allAgents[self.numIDU:(self.numIDU + self.numNIDU)])
        self.ND_agents = deepcopy(allAgents[(self.numIDU + self.numNIDU):])
        self.MSM_agents = []
        self.HIV_agents = []
        self.AIDS_agents = []

        self._set_population()

    def get_agent_charactersitic(self, agent, property):
        """
        :Purpose:
            Determine an agents interal characteristics.
        :Input:	
            agent : int \n
            property: string
        	Allowed properties: \n
        	`HIV`, `Sex Type`, `Drug Type`
        :Output:
            haracteristic : string
        """
        if property in ['HIV', 'Sex Type', 'Drug Type', 'AIDS']:
            try:
                return self.Agents[agent][property]
            except KeyError:
                print self.Agents[agent]
                raise ValueError("Check keys/Agents in self.Agents")
        else:
            raise ValueError("Check agents charactertic! Only 'HIV', 'Sex Type' \
			     'Drug Type' , 'AIDS' possible.")

    def _set_population(self):
        """ 
        :Purpose:
        	Create the population by filling up :py:attr:`self.Agents`.
        	Each agent is a key to an associated dictionary which stores the internal 
        	characteristics in form of an additinoal dictionary of the form
        	``characteristic:value``.
        """
        # Nested dictionary for probability look-up
        # Updated: 1 April 2011
        IDU = {'MSM':{}, 'HM':{}, 'WSW':{}, 'HF':{}}
        IDU['MSM'] = {'HIV':0.67 , 'AIDS':0.058}
        IDU['HM'] = {'HIV':0.40 , 'AIDS':0.058}
        IDU['WSW'] = {'HIV':0.53 , 'AIDS':0.058}
        IDU['HF'] = {'HIV':0.39 , 'AIDS':0.058}
        NIDU = {'MSM':{}, 'HM':{}, 'WSW':{}, 'HF':{}}
        NIDU['MSM'] = {'HIV':0.22 , 'AIDS':0.02}
        NIDU['HM'] = {'HIV':0.048 , 'AIDS':0.002}
        NIDU['WSW'] = {'HIV':0.048 , 'AIDS':0.002}
        NIDU['HF'] = {'HIV':0.048 , 'AIDS':0.002}
        ND = {'MSM':{}, 'HM':{}, 'WSW':{}, 'HF':{}}
        ND['Type'] = ([0,1,2,3], [0.456, 0.485, 0.022, 0.037])
        ND['MSM'] = {'HIV':0.17 , 'AIDS':0.02}
        ND['HM'] = {'HIV':0.015 , 'AIDS':0.0003}
        ND['WSW'] = {'HIV':0.012 , 'AIDS':0.0003}
        ND['HF'] = {'HIV':0.012 , 'AIDS':0.0003}
        ProbLookUp = {'IDU':IDU, 'NIDU':NIDU, 'ND':ND}

        for agent in self.IDU_agents:
            # Sex type
            tmp_rnd = random.random()
            if tmp_rnd < 0.63: SexType = 'HM'
            elif tmp_rnd < 0.63 + 0.249: SexType = 'HF'
            elif tmp_rnd < 0.63 + 0.249 + 0.07: 
                SexType = 'MSM'
                self.MSM_agents.append(agent)
            else: 
                SexType = 'WSW'	
            # HIV  
            prob_HIV = ProbLookUp['IDU'][SexType]['HIV']
            if random.random() < prob_HIV:	
                HIVStatus = 1
                self.HIV_agents.append(agent)
            else:
                HIVStatus = 0
            # AIDS
            prob_AIDS = ProbLookUp['IDU'][SexType]['AIDS']
            if random.random() < prob_AIDS: 
                AIDSStatus = 1
                self.AIDS_agents.append(agent)
            else:
                AIDSStatus = 0
            agent_dict = {'Drug Type': 'IDU'}
            agent_dict.update({'Sex Type':SexType, 'HIV':HIVStatus, 'AIDS':AIDSStatus})
            self.Agents.update({agent:agent_dict})

        for agent in self.NIDU_agents:
            # Sex type
            tmp_rnd = random.random()
            if tmp_rnd < 0.54: SexType = 'HM'
            elif tmp_rnd < 0.54 + 0.332: SexType = 'HF'
            elif tmp_rnd < 0.54 + 0.332 + 0.068: 
                SexType = 'WSW'
            else: 
                SexType = 'MSM'
                self.MSM_agents.append(agent)
            # HIV
            prob_HIV = ProbLookUp['NIDU'][SexType]['HIV']
            if random.random() < prob_HIV:	
                HIVStatus = 1
                self.HIV_agents.append(agent)
            else: 
                HIVStatus = 0
            # AIDS
            prob_AIDS = ProbLookUp['NIDU'][SexType]['AIDS']
            if random.random() < prob_AIDS: 
                AIDSStatus = 1
                self.AIDS_agents.append(agent)
            else:
                AIDSStatus = 0
            agent_dict = {'Drug Type': 'NIDU'}
            agent_dict.update({'Sex Type':SexType, 'HIV':HIVStatus, 'AIDS':AIDSStatus})
            self.Agents.update({agent:agent_dict})

        for agent in self.ND_agents:
            # Sex type
            tmp_rnd = random.random()
            if tmp_rnd < 0.506: SexType = 'HF'
            elif tmp_rnd < 0.506 + 0.453: SexType = 'HM'
            elif tmp_rnd < 0.506 + 0.453 + 0.024: 
                SexType = 'MSM'
                self.MSM_agents.append(agent)
            else:
                SexType = 'WSW'
            # HIV 
            prob_HIV = ProbLookUp['ND'][SexType]['HIV']
            if random.random() < prob_HIV:	
                HIVStatus = 1
                self.HIV_agents.append(agent)
            else:
                HIVStatus = 0
            # AIDS
            prob_AIDS = ProbLookUp['ND'][SexType]['AIDS']
            if random.random() < prob_AIDS: 
                AIDSStatus = 1
                self.AIDS_agents.append(agent)
            else:
                AIDSStatus = 0
            agent_dict = {'Drug Type': 'ND'}
            agent_dict.update({'Sex Type':SexType, 'HIV':HIVStatus, 'AIDS':AIDSStatus})
            self.Agents.update({agent:agent_dict})

    def get_agents(self):
        """ 
        :Purpose:
            Return all agents and their characteristics.

        :Output: 
            :py:attr:`Agents`: dict
        	Dictionary of agents and their characteristics. The agents
        	are the `keys` and a dictionary of `characteristic:value`
        	pair is the entry. 
        """
        return self.Agents

    def print_info(self):
        """ 
        :Purpose:
            Simple fprintf test on std out.
        """
        print ('Number of agents '+str(self.PopulationSize))	
        print ('Number of IDU '+str(len(self.IDU_agents)))	
        print ('Number of NIDU '+str(len(self.NIDU_agents)))	
        print ('Number of ND '+str(len(self.ND_agents)))	
        """
	# random IDU
	agent = random.choice(self.IDU.keys())
	print (' Random IDU agent: ' + str(agent))
	print self.IDU[agent]
	# random NIDU
	agent = random.choice(self.NIDU.keys())
	print (' Random NIDU agent: ' + str(agent))
	print self.NIDU[agent]
	# random ND
	agent = random.choice(self.ND.keys())
	print (' Random ND agent: ' + str(agent))
	print self.ND[agent]
	"""

    def get_info_DrugUserType(self):
        """ 
        :Purpose:
            Return number of IDU, NIDU, and ND users in one array.

        :Output:
            data : dict
        """
        return {'Number of IDU agents': len(self.IDU_agents), 
                'Number of NIDU agents':len(self.NIDU_agents), 
                'Number of ND agents':len(self.ND_agents)}

    def get_info_HIV_IDU(self):
        """
        :Purpose: 
            Return number of HIV among IDU agents. 
            Distinguish between MSM, WSW, HM, and HIF agents.

        :Ooutput: 
            data : dict
        """
        count_HIV_MSM = 0
        count_HIV_HM  = 0
        count_HIV_WSW = 0
        count_HIV_HF  = 0		
        for agent in self.IDU_agents:
            HIVstatus = self.get_agent_charactersitic(agent, 'HIV')
            if HIVstatus == 1:
                SexType = self.get_agent_charactersitic(agent, 'Sex Type')
                if SexType == 'MSM': count_HIV_MSM += 1
                elif SexType == 'HM':count_HIV_HM += 1
                elif SexType == 'WSW':count_HIV_WSW += 1
                elif SexType == 'HF':count_HIV_HF += 1
                else:raise ValueError("Agent must be either HM, MSM, HF, WSW !")
            elif HIVstatus != 0:
                print HIVstatus
                raise ValueError("HIV status must be either 0 or 1 !")
        return {'Number of IDU HIV MSM':count_HIV_MSM, 
                'Number of IDU HIV HM':count_HIV_HM, 
                'Number of IDU HIV WSW':count_HIV_WSW, 
                'Number of IDU HIV HF':count_HIV_HF}

    def get_info_DrugSexType(self):
        """
        :Purpose: 
            Assess the Drug and Sex type prevalences of the population.

        :Ooutput: 
            data : dict
        """
        count_HM	= 0
        count_HF	= 0
        count_MSM	= 0
        count_WSW	= 0
        count_HM_IDU	= 0
        count_HF_IDU	= 0
        count_MSM_IDU	= 0
        count_WSW_IDU	= 0		
        count_HM_NIDU	= 0
        count_HF_NIDU	= 0
        count_MSM_NIDU	= 0
        count_WSW_NIDU	= 0
        count_HM_ND	= 0 
        count_HF_ND	= 0
        count_MSM_ND	= 0
        count_WSW_ND	= 0

        for agent in self.Agents.keys():
            SexType = self.get_agent_charactersitic(agent, 'Sex Type')
            DrugType = self.get_agent_charactersitic(agent, 'Drug Type')
            if SexType == 'HM':
                count_HM += 1 
                if DrugType == 'IDU':count_HM_IDU += 1
                elif DrugType == 'NIDU':count_HM_NIDU += 1
                elif DrugType == 'ND':count_HM_ND += 1
                else: raise ValueError("Drug test must either be IDU, NIDU or ND !")
            elif SexType == 'HF':
                count_HF += 1 
                if DrugType == 'IDU': count_HF_IDU += 1
                elif DrugType == 'NIDU':count_HF_NIDU += 1
                elif DrugType == 'ND':count_HF_ND += 1
                else:raise ValueError("Drug test must either be IDU, NIDU or ND !")
            elif SexType == 'MSM':
                count_MSM += 1 
                if DrugType == 'IDU':count_MSM_IDU += 1
                elif DrugType == 'NIDU':count_MSM_NIDU += 1
                elif DrugType == 'ND':count_MSM_ND += 1
                else:raise ValueError("Drug test must either be IDU, NIDU or ND !")
            elif SexType == 'WSW':
                count_WSW += 1 
                if DrugType == 'IDU':count_WSW_IDU += 1
                elif DrugType == 'NIDU':count_WSW_NIDU += 1
                elif DrugType == 'ND':count_WSW_ND += 1
                else:raise ValueError("Drug test must either be IDU, NIDU or ND !")
            else:raise ValueError("Agent must be either HM, MSM, HF, WSW !")

        return {'Number of HM':count_HM, 
                'Number of HF':count_HF, 
                'Number of MSM':count_MSM, 
                'Number of WSW':count_WSW, 
                'Number of MSM IDU':count_MSM_IDU, 
                'Number of MSM NIDU':count_MSM_NIDU, 
                'Number of MSM ND':count_MSM_ND, 
                'Number of HM IDU':count_HM_IDU, 
                'Number of HM NIDU':count_HM_NIDU, 
                'Number of HM ND':count_HM_ND, 
                'Number of WSW IDU':count_WSW_IDU, 
                'Number of WSW NIDU':count_WSW_NIDU, 
                'Number of WSW ND':count_WSW_ND,
                'Number of HF IDU':count_HF_IDU, 
                'Number of HF NIDU':count_HF_NIDU, 
                'Number of HF ND':count_HF_ND}

class TestClassMethods(unittest.TestCase):
    """ 
    :Purpose:
        unittest
    """
    def setUp(self):
        """
        :Purpose:
            Tests that all models from setup pass inspection. ``setUp`` is perfomed before each method.
        """
        self.N_pop = 10000

    def test_HIV(self):
        """ Test: Testing HIV agent array"""
        print "\t__Testing the HIV agent list"	
        tmpCount = 0
        myPopulation = PopulationClass(n=self.N_pop)
        for a in myPopulation.Agents.keys():
            HIVstatus = myPopulation.get_agent_charactersitic(a,'HIV')
            if HIVstatus != 0:
                tmpCount += 1
                self.assertTrue(a in myPopulation.HIV_agents)
        self.assertTrue(len(myPopulation.HIV_agents) == tmpCount)

    def test_AIDS(self):
        """ Test: Testing AIDS agent array"""
        print "\t__Testing the AIDS agent list"	
        tmpCount = 0
        myPopulation = PopulationClass(n=self.N_pop)
        for a in myPopulation.Agents.keys():
            AIDSstatus = myPopulation.get_agent_charactersitic(a,'AIDS')
            if AIDSstatus != 0:
                tmpCount += 1
                self.assertTrue(a in myPopulation.AIDS_agents)
        self.assertTrue(len(myPopulation.AIDS_agents) == tmpCount)

    def test_consistency(self):
        """ Test: Testing consistency"""
        print "\t__Testing consistency of agent lists"
        myPop = PopulationClass(n=self.N_pop)
        NormalAgents = list(set(range(myPop.PopulationSize)).difference(
            set(myPop.IDU_agents).union(set(myPop.MSM_agents))))
        MSMAgents = list(set(myPop.MSM_agents).difference(set(myPop.IDU_agents)))
        IDUagents = myPop.IDU_agents
        NumAgents = len(NormalAgents)+len(MSMAgents)+len(IDUagents)
        self.assertTrue(NumAgents == self.N_pop,"NumAgents=%d, \
				PopulationSize = %d"%(NumAgents, self.N_pop))

    def test_MSM(self):
        """ Test: Testing MSM agent array"""
        print "\t__Testing the MSM agent list"	
        tmpCount = 0
        myPopulation = PopulationClass(n=self.N_pop)
        for a in myPopulation.Agents.keys():
            SEXstatus = myPopulation.get_agent_charactersitic(a,'Sex Type')
            if SEXstatus == 'MSM':
                tmpCount += 1
                self.assertTrue(a in myPopulation.MSM_agents)
        self.assertTrue(len(myPopulation.MSM_agents) == tmpCount)

    def test_Obesity(self):
        """ Test: Testing the population"""
        print "\t__Testing the population"	
        myPopulation = PopulationClass(n=self.N_pop)
        for agent in myPopulation.Agents.keys():
            tmp_DrugType = myPopulation.get_agent_charactersitic(agent,'Drug Type')
            self.assertTrue(tmp_DrugType in ['NIDU','IDU','ND'])
            if tmp_DrugType == 'NIDU':
                self.assertTrue(agent in myPopulation.NIDU_agents)
            elif tmp_DrugType == 'IDU':
                self.assertTrue(agent in myPopulation.IDU_agents)
            else: self.assertTrue(agent in myPopulation.ND_agents)
        self.assertTrue(len(myPopulation.NIDU_agents) == myPopulation.numNIDU)
        self.assertTrue(len(myPopulation.IDU_agents) == myPopulation.numIDU)
        self.assertTrue(len(myPopulation.ND_agents) == myPopulation.numND)

if __name__=='__main__':
    """unittest"""
    unittest.main()				