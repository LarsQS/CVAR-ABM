#!/usr/bin/env python
# encoding: utf-8
"""
Module that carries the defintions for the Population for 
the ABM of drug use and HIV.

+ Author:		Lars Seemann
+ Email:		lseemann@uh.edu
+ Date:		2011-01-28

Copyright (c) 2010, under the Simplified BSD License. \n
For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php \n
All rights reserved.
"""
__author__="Lars Seemann (lseemann@uh.edu)"

import random
import unittest
from copy import deepcopy, copy

try: from HIVABM_PartialNetwork import ScaleFreeNetworkClass
except ImportError: raise ImportError("Can't import ScaleFreeNetworkClass")


class HIVModel(ScaleFreeNetworkClass):
    """
    :Purpose:
    	This is the core class used to simulate 
    	the spread of HIV and drug use in one MSA 
    	(Metropolitan Statistical Area).

    :Input:	
    	N : int
    		Number of agents. Default: 10000

    	tmax: int	
    		Number of simulation steps (years).

    	:py:class:`SocialNetworkClass` : Inherited

    	:py:class:`PopulationClass` : Inherited

    :Attributes:
    	:py:attr:`tmax` : int
    		Number of time steps simulated.

    	:py:attr:`CleanSyringeUsers` : list

    	:py:attr:`SEPAgents` : dict
            Dictionary of users who participated in a 
            syringe exchange program (SEP) {agent:time}.

    	:py:attr:`DrugTreatmentAgents` : dict
            Dictionary of users who underwent drug 
            treatment {agent:time}.

    	:py:attr:`TestedAgents` : list
            List of agents who get tested for HIV every time step.

    	:py:attr:`tmp_Agents` : dict
            Changes resulting from parsing through the agents 
            and applying the update rules are stored 
            in :py:attr:`tmp_agent_dict`.

    	All attributes from :py:class:`SocialNetworkClass` \n

    	All attributes from :py:class:`PopulationClass` 

    :Methods:
    	:py:meth:`_update_population` \n
    	:py:meth:`_needle_transmission` \n
    	:py:meth:`_sex_transmission` \n
    	:py:meth:`_drug_transition` \n
    	:py:meth:`_update_IDU` \n
    	:py:meth:`_update_NIDU_ND` \n
    	:py:meth:`_update_AllAgents` \n
    	:py:meth:`_VCT` \n
        :py:meth:`_SEP` \n
        :py:meth:`_enter_drug_treatment` \n
        :py:meth:`_initiate_HAART` \n
        :py:meth:`_get_partner` \n
        :py:meth:`_update_AllAgents` \n
    	:py:meth:`run` \n
    	:py:meth:`store_results` \n
    	:py:meth:`get_HIV_prevalence_drugs` \n
    	:py:meth:`get_HIV_prevalence` \n
    	:py:meth:`plot_results` \n
    	All methods from :py:class:`SocialNetworkClass` \n
    	All methods from :py:class:`PopulationClass`
    """

    def __init__(self, N = 10000, tmax = 40):
        """ Initialize HIVModel object """
        if (type(tmax) is not int):
            raise ValueError("Number of time steps must be integer")
        else:
            self.tmax = tmax
        self.N_TreatmentSpots = 0    # initiate total number of treatment spots

        # Initialize inherited social Network
        #print 'Create Social Network...'
        #PopulationClass.__init__(self, n = N)	# Create population

        # OLD
        # Clean syring users
        #self.CleanSyringeUsers =[]
        #for agent in self.Agents.keys():
        #    agent_drug_type = self.get_agent_charactersitic(agent, 'Drug Type')
        #    if agent_drug_type == 'IDU' and random.random() < 0.1:
        #        self.CleanSyringeUsers.append(agent)

        # Other lists / dictionaries
        self.SafeSexUsers = dict()			# dict of users who had safe sex (agent:time)
        self.SEPAgents = dict()				# dict of users who used SEP ( agent:time)
        self.DrugTreatmentAgents_current = dict()	# dictionary of users who are currently undergoing  
                                                        # drug treatent (agent:time)
        self.DrugTreatmentAgents_past = dict()	        # dictionary of users who underwent drug treat- 
                                                        # ment in the past(agent:time)
        self.VCTAgents = dict()			        # list of agents who get tested for HIV ((agent:time) 
        self.HAARTAgents = dict()			# dict of agents in HIV treatment (agent:time) 
        self.HAARTAdherence = dict()			# dict of agents HIV treatment (HAART) adherence (agent:adherence)
        self.tmp_Agents = dict()			# dictionary that keeps track of the update 


        # Keep track of results
        self.IDUprevalence = [0]*self.tmax			
        self.NIDUprevalence = [0]*self.tmax			
        self.NDprevalence = [0]*self.tmax
        self.HIVtotalPrevalence = [0]*self.tmax				

        self.HIVinIDU = [0]*self.tmax	
        self.HIVinNIDU = [0]*self.tmax	
        self.HIVinND = [0]*self.tmax	

    def _update_population(self):
        """
        :Purpose:
        	Update the population. Changes resulting from parsing
        	through the agents and applying the update rules are stored 
        	in :py:attr:`tmp_agent_dict`. This method updates the whole 
        	population, i.e., it copies changes from the :py:attr:`tmp_agent_dict`
        	dictionary and copies it into the :py:attr:`Agents` dictionary.
        :Input:
            none

        :Output: 
            none
        """
        for (u,d) in self.tmp_Agents.iteritems():
            self.Agents.update({u:d})

    def _needle_transmission(self, agent, partner):
        """ 
        :Purpose:
        	Simulate random transmission of HIV between two IDU agents through needle.\n 
        	Needed in _update_IDUand
        :Input:
        	agents : int
        	partner : int
        :Output: -
        """
        # both must be IDU
        partner_drug_type = self.get_agent_charactersitic(partner, 'Drug Type')
        agent_drug_type = self.get_agent_charactersitic(agent, 'Drug Type')
        if not (partner_drug_type == 'IDU' and agent_drug_type == 'IDU' ):
            raise ValueError("To share a needle both agents must be IDU")
        else:
            # Do they share a needle?
            if (agent in self.SEPAgents.keys() or partner in self.SEPAgents.keys() ):
                pass # no needle sharing
            elif (random.random() < 0.66): 
                pass # no needle sharing
            else:								# they do share a needle
                # HIV+ ?
                HIV_agent = self.get_agent_charactersitic(agent,'HIV')
                HIV_partner = self.get_agent_charactersitic(partner,'HIV')
                # if agent HIV+
                if HIV_agent == 1 and random.random() < 0.009:
                    tmp_partner_dict = deepcopy(self.Agents[partner])
                    tmp_partner_dict.update({'HIV' : 1})
                    self.tmp_Agents.update({partner : tmp_partner_dict})
                # if partner HIV+
                if HIV_partner == 1 and random.random() < 0.009:
                    tmp_agent_dict = deepcopy(self.Agents[agent])
                    tmp_agent_dict.update({'HIV' : 1})
                    self.tmp_Agents.update({agent : tmp_partner_dict})

    def _sex_transmission(self, agent, partner, time):
        """ 
        :Purpose:
            Simulate random transmission of HIV between two agents through Sex. 
            Needed for all users. Sex is not possible in case the agent and 
            assigned partner have incompatible Sex behavior.

        :Input:
            agents : int
            partner : int
            time : int

        :Output:
            none
        """
        
        # Double check: Sex possible?
        Type_agent   = self.get_agent_charactersitic(agent,'Sex Type')
        Type_partner = self.get_agent_charactersitic(partner,'Sex Type')
        if Type_agent == 'HM' and Type_partner in [ 'HF', 'MSM']:
            pass
        elif Type_partner == 'HM' and Type_agent in [ 'HF', 'MSM']:
            pass
        elif Type_agent == 'MSM' and Type_partner in ['MSM', 'WSW' , 'HF']:
            pass
        elif Type_partner == 'MSM' and Type_agent in ['MSM', 'WSW' , 'HF']:
            pass
        elif Type_agent == 'WSW' and Type_partner in ['MSM', 'WSW' , 'HM']:
            pass
        elif Type_partner == 'WSW' and Type_agent in ['MSM', 'WSW' , 'HM']:
            pass
        else:
            raise ValueError("Sex must be possible! %s %s"%(str(Type_agent),str(Type_partner)))
        
        # HIV status of agent and partner
        # Everything from here is only run if one of them is HIV+
        HIVstatus_Agent = self.get_agent_charactersitic(agent,'HIV')
        HIVstatus_Partner = self.get_agent_charactersitic(partner,'HIV')
        
        if HIVstatus_Agent == 1 or HIVstatus_Partner == 1:
            # Sex between men?
            if Type_agent in ['HM', 'MSM'] and Type_partner in ['HM', 'MSM']:
                SexBetweenMen = 1
            elif Type_partner in ['HM', 'MSM'] and Type_agent in ['HM', 'MSM']:
                SexBetweenMen = 1
            else:
                SexBetweenMen = 0
            # Define probabilities for unsafe sex
            UnsafeSexProb_UnknowHIVStatus = {'MSM':{'IDU':0.65, 'NIDU':0.55, 'ND':0.45},
                                             'HM':{'IDU':0.75, 'NIDU':0.75, 'ND':0.75},
                                             'HF':{'IDU':0.75, 'NIDU':0.75, 'ND':0.75},
                                             'WSW':{'IDU':0.80, 'NIDU':0.80, 'ND':0.80}}
            UnsafeSexProb_HIVPosStatus = {'MSM':{'IDU':0.65, 'NIDU':0.55, 'ND':0.45},
                                             'HM':{'IDU':0.40, 'NIDU':0.40, 'ND':0.40},
                                             'HF':{'IDU':0.40, 'NIDU':0.40, 'ND':0.40},
                                             'WSW':{'IDU':0.45, 'NIDU':0.45, 'ND':0.45}}
            DrugType_Agent = self.get_agent_charactersitic(agent,'Drug Type')
            DrugType_Partner = self.get_agent_charactersitic(partner,'Drug Type')
            if agent in self.VCTAgents and HIVstatus_Agent == 1:
                p_SafeSex_Agent = UnsafeSexProb_HIVPosStatus[Type_agent][DrugType_Agent]
            else:
                p_SafeSex_Agent = UnsafeSexProb_HIVPosStatus[Type_agent][DrugType_Agent]
            if partner in self.VCTAgents and HIVstatus_Partner == 1:
                p_SafeSex_Partner = UnsafeSexProb_HIVPosStatus[Type_partner][DrugType_Partner]
            else:
                p_SafeSex_Partner= UnsafeSexProb_HIVPosStatus[Type_partner][DrugType_Partner]
            p_SafeSex = 0.5*(p_SafeSex_Agent + p_SafeSex_Partner)
            # Safe Sex?
            if random.random() < p_SafeSex:
                SafeSexFlag = 0
            else:
                SafeSexFlag = 1
            # HIV transmission only possible if unprotected sex
            if SafeSexFlag == 0:
                # HIV transmission probability
                if agent in self.HAARTAgents:
                    if SexBetweenMen == 1:
                        if random.random() < 0.4:
                            p_transmission = random.choice([0.452, 0.382, 0.214, 0.113])
                        else:
                            p_transmission = 0.010
                    elif SexBetweenMen == 0:
                        if random.random() < 0.4:
                            p_transmission = random.choice([0.113, 0.092, 0.047, 0.024, 0.002])
                        else:
                            p_transmission = 0.002
                    else:
                        raise ValueError("Check SexBetweenMen flag! %s"%str(SexBetweenMen))
                else:
                    if SexBetweenMen == 1:
                        p_transmission = 0.452
                    elif SexBetweenMen == 0:
                        p_transmission = 0.113
                    else:
                        raise ValueError("Check SexBetweenMen flag! %s"%str(SexBetweenMen))
                # HIV transmission
                # if agent HIV+
                if HIVstatus_Agent == 1 and random.random() < p_transmission:
                    tmp_partner_dict = deepcopy(self.Agents[partner])
                    tmp_partner_dict.update({'HIV' : 1})
                    self.tmp_Agents.update({partner : tmp_partner_dict})
                # if partner HIV+
                if HIVstatus_Partner == 1 and random.random() < p_transmission:
                    tmp_agent_dict = deepcopy(self.Agents[agent])
                    tmp_agent_dict.update({'HIV' : 1})
                    self.tmp_Agents.update({agent : tmp_agent_dict})

    def _drug_transition(self, agent, partner):
        """ 
        :Purpose:
            Simulate transition of drug behavior. The following scenarios are 
            possible:
            + ND agent might become NIDU when meeting NIDU 
            + NIDU might become IDU when meeting IDU 
            The function is only applied for NIDU and ND users.

        :Input:
            agents : int
            partner : int

        :Output: -
        """
        partner_drug_type = self.get_agent_charactersitic(partner, 'Drug Type')
        agent_drug_type = self.get_agent_charactersitic(agent, 'Drug Type')
        # NIDU -> IDU
        if agent_drug_type == 'NIDU' and partner_drug_type == 'IDU':
            if random.random() < 0.12:
                # agent becomes IDU
                tmp_agent_dict = deepcopy(self.Agents[agent])
                tmp_agent_dict.update({'Drug Type' : 'IDU'})
                self.tmp_Agents.update({agent : tmp_agent_dict})
        if  partner_drug_type == 'NIDU' and agent_drug_type == 'IDU':
            if random.random() < 0.12:
                # partner becomes IDU
                tmp_partner_dict = deepcopy(self.Agents[partner])
                tmp_partner_dict.update({'Drug Type' : 'IDU'})
                self.tmp_Agents.update({partner : tmp_partner_dict})	
        # ND -> NIDU
        if agent_drug_type == 'ND' and partner_drug_type == 'NIDU':
            if random.random() < 0.30:
                # agent becomes NIDU
                tmp_agent_dict = deepcopy(self.Agents[agent])
                tmp_agent_dict.update({'Drug Type' : 'NIDU'})
                self.tmp_Agents.update({agent : tmp_agent_dict})
        if partner_drug_type == 'ND' and agent_drug_type == 'NIDU':
            if random.random() < 0.30:
                # partner becomes NIDU
                tmp_partner_dict = deepcopy(self.Agents[partner])
                tmp_partner_dict.update({'Drug Type' : 'NIDU'})
                self.tmp_Agents.update({partner : tmp_partner_dict})	

    def _update_IDU(self, agent, partner, time):
        """
        :Purpose:
            Let IDU agent interact with a partner.
            Update IDU agents:
                    1 - determine transition type
                    2 - Injection rules
                    3 - Sex rules
                    4 - HIV transmission
                    5 - SEP

        :Input:	
            agent : int

            partner : int 

            time : int

        Output:	
            none

        """
        partner_drug_type = self.get_agent_charactersitic(partner, 'Drug Type')
        agent_drug_type = self.get_agent_charactersitic(agent, 'Drug Type')
        if partner_drug_type == 'IDU' and agent_drug_type == 'IDU':
            if random.random() < 0.5: self. _needle_transmission(agent, partner)
            else: self._sex_transmission(agent, partner, time)
        elif partner_drug_type in ['NIDU', 'ND'] or agent_drug_type in ['NIDU', 'ND']:
            self._drug_transition(agent, partner)
            self._sex_transmission(agent, partner, time)
        else:
            raise ValueError("Agents must be either IDU, NIDU, or ND")

        # SEP
        if partner_drug_type == 'IDU':
            if random.random() < 0.6: self.SEPAgents.update({partner:time})
        if agent_drug_type == 'IDU':
            if random.random() < 0.6: self.SEPAgents.update({agent:time})

    def _update_NIDU_ND(self, agent, partner, time):
        """
        :Purpose:
            Let NIDU or ND agent interact.

        :Input:	

            agent : int
            partner : int 
            time : int

        Output: none
        """	
        self._drug_transition(agent, partner)
        self._sex_transmission(agent, partner, time)

    def _VCT(self, agent, time):
        """
        :Purpose:
            Account for voluntary Counseling and Testing(VCT)

        :Input:	
            agent : int
            partner : int
            time : int

        :Output:
            none
        """	
        # Drug Treatment
        # SEP
        drug_type = self.get_agent_charactersitic(agent, 'Drug Type')
        SEPstat = FALSE
        if agent in self.SEPAgents:
            if time == self.SEPAgents[agent]:
                SEPstat = TRUE
            # MSM
        if drug_type == 'IDU':
            if SEPstat:
                if random.random() < 0.75:
                    self.VCTAgents.update({agent:time})
            else:
                if random.random() < 0.60:
                    self.VCTAgents.update({agent:time})
        elif agent in self.MSM_agents and random.random() < 0.25:
            self.VCTAgents.update({agent:time})
        else:
            if random.random() < 0.09:
                self.VCTAgents.update({agent:time})

    def _SEP(self, agent, time):
        """
        :Purpose:
            Account for SEP (Syringe Exchange Program) for IDU agents. \n
            SEP use depends on the year and the functional relationship
            is given by
            P(SEP USE (t+1) | IDU) = P(SEP USE (t) | IDU)*f_t(SEP density)
            where f_t(SEP density) = (1+? + ln(?))

        :Input:	
            time : int

        :Output:
            bool
        """
        if agent not in self.IDU_agents:
            raise ValueError("_SEP only valid for IDU agents! %s"%
                             str(self.get_agent_charactersitic(agent,'Drug Type')))
        # The following values are estimated from the SEP prob plot in the documentation
        P_SEPuse_dict={1:0.25, 2:0.275, 3:0.3, 4:0.32, 5:0.34, 6:0.36, 7:0.38, 8:0.39, 9:0.41, 10:0.42, 11:0.44}
        if random.random()<P_SEPuse_dict[time]:
            self.SEPAgents.update({agent:time})

    def _enter_drug_treatment(self, agent, time):
        """
        :Purpose:
            Account for drug treatment. \n
            Entering drug treatment is similar to SEP, drug treatment has a functional relationship given as 
            follows:
            P(treatment (t+1) | IDU or NIDU) = P(treatment (agent,t) | IDU) if x < N
            else: 0 (x >= 0)
            where N is the total number of treatment slots available. \n
            An agent who was already in drug treatment and relapsed, has a pobability twice as strong
            to reenter drug treatment at a later point.

        :Input:	
            agent : int

            time : int

        :Output:
            bool
        """
        N_TrSpots_Max = 500   # max number of treatment spots
        CurrentlyTreatedFlag = FALSE

        if self.N_TreatmentSpots < 500:
            # Determine probability of treatment
            if agent in self.DrugTreatmentAgents_current.keys():
                CurrentlyTreatedFlag = TRUE
                prob = 0.6
            elif agent in self.DrugTreatmentAgents_past.keys():
                if agent in self.SEPAgents:
                    prob = 2.*0.09
                elif agent in self.IDU_agents or agent in self.NIDU_agents:
                    prob = 2.*0.18
            else:
                # Probability considering SEP status
                if agent in self.SEPAgents:
                    prob = 0.09
                elif agent in self.IDU_agents or agent in self.NIDU_agents:
                    prob = 0.18
                else:
                    raise ValueError("_SEP only valid for NIDU or IDU agents! %s"%
                                     str(self.get_agent_charactersitic(agent,'Drug Type'))) 
            # Assign treatment randomly
            if random.random() < factor*0.09:
                self.DrugTreatmentAgents_current.update({agent:time})
                self.N_TreatmentSpots += 1
                return TRUE
            else:
                if CurrentlyTreatedFlag:
                    # Exit drug treatment
                    self.DrugTreatmentAgents_past.update({agent:self.DrugTreatmentAgents_current[agent]})
                    del self.DrugTreatmentAgents_current[agent]
                return FALSE
        elif self.N_TreatmentSpots == 500:
            return FALSE
        else:
            raise ValueError("Check self.N_TreatmentSpots! Max value = 500! %d"%self.N_TreatmentSpots)

    def _initiate_HAART(self, agent, time):
        """
        :Purpose:
            Account for HIV treatment through highly active antiretroviral therapy (HAART).
            HAART was implemented in 1996, hence, there is treatment only after 1996.
            HIV treatment assumes that the agent knows his HIV+ status.

        :Input:	
            time : int

        :Output:
            none
        """
        # Determine probability of HIV treatment
        if time >= 5:    # HAART implemented 1996
            if agent in self.HAARTAgents:
                pass    # agents never discontinure therapy
            if agent in self.VCTAgents:
                agent_drug_type = self.get_agent_charactersitic(agent, 'Drug Type')
                if agent_drug_type == 'IDU':
                    if agent in self.VCTAgents:
                        if agent in self.DrugTreatmentAgents_current:
                            prob = 0.21
                        else:
                            prob = 0.12
                elif agent_drug_type == 'NIDU':
                    if agent in self.VCTAgents:
                        prob = 0.14
                else:
                    if agent in self.VCTAgents:
                        prob = 0.22
            else:
                prob = 0.0
        else:
            prob = 0.0

        if random.random() < prob:
            self.HAARTAgents.update({agent:time})
            # adherence
            if random.random() < 0.6:
                adherence = 1
            else:
                adherence = 0
            self.HAARTAdherence.update({agent:adherence})

    def _get_partner(self, agent, agent_drug_type):
        """
        :Purpose:
            Get a partner for agent.

        :Input:	
            agent : int

            agent_drug_type : str
              Either 'IDU', 'NIDU', 'ND'

        :Output:
            partner : int
        """
        if agent not in range(self.PopulationSize):
            raise ValueError("Invalid agent! %s"%str(agent))
        if agent_drug_type not in ['IDU', 'NIDU', 'ND']:
            raise ValueError("Invalid drug type! %s"%str(agent_drug_type))

        if agent_drug_type == 'NIDU':
            pass
        elif agent_drug_type == 'IDU':
            pass
        else:
            pass
        
    def _get_number_of_partners(self,
                                agent,
                                agent_drug_type,
                                agent_sex_type):
        """
        :Purpose:
            Get number of partners for a agent. 
            Drawn from Poisson distribution.

        :Input:	
            agent : int

            agent_drug_type : str
              Either 'IDU', 'NIDU', 'ND'

        :Output:
            NumPartners : int
        """
        # Check input
        # Drug type
        if agent not in range(self.PopulationSize):
            raise ValueError("Invalid agent! %s"%str(agent))
        if agent_drug_type not in ['IDU', 'NIDU', 'ND']:
            raise ValueError("Invalid drug type! %s"%str(agent_drug_type))
        # Sex type
        if agent_sex_type not in ['HM','NF','MSM','WSW']:
            raise ValueError("Invalid sex type! %s"%str(agent_sextype))
        if agent_drug_type == 'NIDU':
            PoissonLambda = 3
        elif agent_drug_type == 'IDU':
            PoissonLambda = 5
        elif agent_sex_type == 'MSM':
            PoissonLambda = 1 
        RandNumCont = np.random.poisson(PoissonLambda, 1)
        return RandNumCont[0]
    
    def _update_AllAgents(self, time):
        """
        :Purpose:
            Update IDU agents:\n
            For each agent:\n
                1 - determine agent type
                2 - update accordingly
                3 - VCT (Voluntsry Counseling and Testing)

        :Input:
            agent, time

        :Output:
            none

        """	
        for agent in self.Agents.keys():
            agent_drug_type = self.get_agent_charactersitic(agent, 'Drug Type')
            NumPartners = self._get_number_of_partners(agent_drug_type)
            n = 0
            while n < NumPartners:
                # Choose a random partner
                partner = self._get_partner(agent, agent_drug_type)
                # Transmission of HIV / Transition of agents
                if agent_drug_type == 'IDU':
                    self._update_IDU(agent, partner, time)
                else:
                    self._update_NIDU_ND(agent, partner, time)
                # Counseling and Testing (VCT)
                self._VCT(agent,time)
                # 

    def run(self):
        """ 
        Core of the model 
        :param agent: agent
        """
        for t in range(1,self.tmax):
            print ('TIME = '+str(t)+' ________________________________________')
            self._update_AllAgents(t)
            self._update_population()
            self.store_results(t)

    def store_results(self, time):
        """ Store result into appropriate arrays """
        NumIDU = 0
        NumNIDU = 0
        NumND = 0
        for agent in self.Agents.keys():
            agent_drug_type = self.get_agent_charactersitic(agent, 'Drug Type')
            if agent_drug_type == 'IDU': NumIDU+=1
            elif agent_drug_type == 'NIDU': NumNIDU+=1
            elif agent_drug_type == 'ND': NumND+=1
            else: raise ValueError("Agents must be either IDU, NIDU, or ND")

        self.IDUprevalence[time] = NumIDU
        self.NIDUprevalence[time] = NumNIDU
        self.NDprevalence[time] = NumND
        self.HIVtotalPrevalence[time] = self.get_HIV_prevalence()

        HIVinDrugsResults = self.get_HIV_prevalence_drugs()
        self.HIVinIDU[time] = HIVinDrugsResults[0]
        self.HIVinNIDU[time] = HIVinDrugsResults[1]
        self.HIVinND[time] = HIVinDrugsResults[2]

    def get_HIV_prevalence_drugs(self):
        """ 
        get HIV prevalence within all three drug user groups
        """
        count_HIV_IDU = 0
        count_HIV_NIDU  = 0
        count_HIV_ND = 0

        for agent in self.Agents:
            HIVstatus = self.get_agent_charactersitic(agent, 'HIV')
            if HIVstatus == 1:
                agent_drug_type = self.get_agent_charactersitic(agent, 'Drug Type')
                if agent_drug_type == 'IDU': count_HIV_IDU += 1
                elif agent_drug_type == 'NIDU': count_HIV_NIDU += 1
                elif agent_drug_type == 'ND': count_HIV_ND += 1
                else: raise ValueError("Agent must be either IDU, NIDU or ND !")
            elif HIVstatus != 0:
                print HIVstatus
                raise ValueError("HIV status must be either 0 or 1 !")
        print [count_HIV_IDU, count_HIV_NIDU, count_HIV_ND]
        return [count_HIV_IDU, count_HIV_NIDU, count_HIV_ND]

    def get_HIV_prevalence(self):
        """ get HIV prevalence within IDUs """
        HIVcount = 0.0
        for agent in self.Agents.keys():
            HIVstatus = self.get_agent_charactersitic(agent, 'HIV')
            if HIVstatus == 1: HIVcount += 1
        print HIVcount
        return HIVcount


    def plot_results(self):
        """ Plot Results """
        plt.subplot(231)
        plt.plot(self.IDUprevalence[1:], 'bo-', markerfacecolor='green')
        plt.ylabel('IDU')
        plt.xlabel('Time')

        plt.subplot(232)
        plt.plot(self.NIDUprevalence[1:], 'bo-', markerfacecolor='green')
        plt.ylabel('NIDU')
        plt.xlabel('Time')

        plt.subplot(233)
        plt.plot(self.NDprevalence[1:], 'bo-', markerfacecolor='green')
        plt.ylabel('ND')
        plt.xlabel('Time')

        plt.subplot(234)
        plt.plot(self.HIVinIDU[1:], 'bo-', markerfacecolor='red')
        plt.ylabel('IDU HIV+')
        plt.xlabel('Time')

        plt.subplot(235)
        plt.plot(self.HIVinNIDU[1:], 'bo-', markerfacecolor='red')
        plt.ylabel('NIDU HIV+ ')
        plt.xlabel('Time')

        plt.subplot(236)
        plt.plot(self.HIVinND[1:], 'bo-', markerfacecolor='red')
        plt.ylabel('ND HIV+ ')
        plt.xlabel('Time')

        plt.show()



class TestClassMethods(unittest.TestCase):
    """ 
    :Purpose:
    	Unittest for testing the methods of the HIV class.
    """
    def setUp(self):
        """
        :Purpose:
            Tests that all models from setup pass inspection. 
            ``setUp`` is perfomed before each method.
        """
        self.N_pop = 1000
        self.T = 40

    def test_update_population(self):
        """ Test: :py:meth:`_update_population` """
        print " ... Test: method _update_population"
        myModel = HIVModel(N = self.N_pop, tmax = self.T)
        self.assertTrue(len(myModel.tmp_Agents) == 0)
        myModel.tmp_Agents = {1:{'N':1, 'M':2}, 2:{'N':3, 'M':4}}
        myModel._update_population()
        self.assertTrue(myModel.Agents[1]['N'] == 1)
        self.assertTrue(myModel.Agents[1]['M'] == 2)
        self.assertTrue(myModel.Agents[2]['N'] == 3)
        self.assertTrue(myModel.Agents[2]['M'] == 4)

    def test_needle_transmission(self):
        """ Test: :py:meth:`_needle_transmission` """
        print " ... Test: method _needle_transmission"	
        myModel = HIVModel(N = self.N_pop, tmax = self.T)
        self.assertTrue(len(myModel.tmp_Agents) == 0)
        chosenpairs = {}
        for agent in myModel.IDU_agents:
            partner = random.choice(myModel.IDU_agents)
            chosenpairs.update({agent:partner})
            myModel._needle_transmission(agent, partner)

        if len(myModel.tmp_Agents) > 0:
            for (agent, d) in myModel.tmp_Agents.iteritems():
                self.assertTrue(agent in myModel.IDU_agents)
                self.assertTrue(d['HIV'] == 1)
                # agent or partner must have been HIV=1 before transmission
                HIVstatus = (myModel.Agents[chosenpairs[agent]]['HIV'] == 1
                             or myModel.Agents[agent]['HIV'] == 1)
                self.assertTrue(HIVstatus == 1)

    def test_sex_transmission(self):
        """ Test: :py:meth:`_sex_transmission` """
        print " ... Test: method _sex_transmission"	
        myModel = HIVModel(N = self.N_pop, tmax = self.T)
        self.assertTrue(len(myModel.tmp_Agents) == 0)
        chosenpairs = {}
        for agent in myModel.IDU_agents:
            partner = random.choice(myModel.IDU_agents)
            chosenpairs.update({agent:partner})
            time = random.choice(range(self.T))
            myModel._sex_transmission(agent, partner, time)

        if len(myModel.tmp_Agents) > 0:
            for (agent, d) in myModel.tmp_Agents.iteritems():
                self.assertTrue(d['HIV'] == 1)
                # agent or partner must have been HIV=1 before transmission
                HIVstatus = (myModel.Agents[chosenpairs[agent]]['HIV'] == 1
                             or myModel.Agents[agent]['HIV'] == 1)
                self.assertTrue(HIVstatus == 1)
                # check if sex behavior matches
                Type_agent = myModel.get_agent_charactersitic(agent, 'Sex Type')
                Type_partner = myModel.get_agent_charactersitic(partner, 'Sex Type')
                Type_agent   = self.get_agent_charactersitic(agent,'Sex Type')
                Type_partner = self.get_agent_charactersitic(partner,'Sex Type')
                if Type_agent == 'HM': self.assertTrue(Type_partner in [ 'HF', 'MSM'])
                elif Type_partner == 'HM':  self.assertTrue(Type_agent in [ 'HF', 'MSM'])
                elif Type_agent == 'MSM': self.assertTrue(Type_partner in ['MSM', 'WSW' , 'HF'])
                elif Type_partner == 'MSM': self.assertTrue(Type_agent in ['MSM', 'WSW' , 'HF'])
                elif Type_agent == 'WSW': self.assertTrue(Type_partner in ['MSM', 'WSW' , 'HM'])
                elif Type_partner == 'WSW': self.assertTrue(Type_agent in ['MSM', 'WSW' , 'HM'])


if __name__=='__main__':
    """unittest"""
    unittest.main()	

