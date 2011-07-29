#!/usr/bin/env python
# encoding: utf-8
"""
HIV_AssessPopulation.py

Purpose:	Assess the population created by HIV_PopulationClass.py

Output:		List of popluation characteristics on std out

Author:		Lars Seemann
Email:		lseemann@uh.edu
Date:		2010-08-19

Copyright (c) 2010, under the Simplified BSD License.
For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php
All rights reserved.
"""

import random
import copy
from numpy import array, zeros
#import matplotlib.pyplot as plt

from HIVABM_Population import PopulationClass

# number of runs
N = 100
PopulationSize = 10000

# arrays for results

# Drug type
IDU = zeros(N)
NIDU = zeros(N)
ND = zeros(N)

# Drug type by sex be= zeros(N)havior
HM = zeros(N)
HF = zeros(N)
MSM = zeros(N)
WSW = zeros(N)

HM_IDU = zeros(N)
HF_IDU = zeros(N)
MSM_IDU = zeros(N)
WSW_IDU = zeros(N)

HM_NIDU = zeros(N)
HF_NIDU = zeros(N)
MSM_NIDU = zeros(N)
WSW_NIDU = zeros(N)

HM_ND = zeros(N)
HF_ND = zeros(N)
MSM_ND = zeros(N)
WSW_ND = zeros(N)

# HIV
# HIV cases within IDUs
IDU_HIV_MSM = zeros(N)
IDU_HIV_HM = zeros(N)
IDU_HIV_WSW = zeros(N)
IDU_HIV_HF = zeros(N)


for i in range(N):
    myModel = PopulationClass(n=PopulationSize)

    # Drug User
    DrugUser = myModel.get_info_DrugUserType()
    IDU[i] = DrugUser['Number of IDU agents'] / float(PopulationSize)
    NIDU[i] = DrugUser['Number of NIDU agents'] / float(PopulationSize)
    ND[i] = DrugUser['Number of ND agents'] / float(PopulationSize)

    # Drug user by sex behavior 
    DrugSexUser = myModel.get_info_DrugSexType()
    HM[i] = DrugSexUser['Number of HM'] / float(PopulationSize)
    HF[i] = DrugSexUser['Number of HF'] / float(PopulationSize)
    MSM[i] = DrugSexUser['Number of MSM'] / float(PopulationSize)
    WSW[i] = DrugSexUser['Number of WSW'] / float(PopulationSize)

    HM_IDU[i] = DrugSexUser['Number of HM IDU'] / float( DrugUser['Number of IDU agents'])
    HF_IDU[i] = DrugSexUser['Number of HF IDU'] / float( DrugUser['Number of IDU agents'])
    MSM_IDU[i] = DrugSexUser['Number of MSM IDU'] / float( DrugUser['Number of IDU agents'])
    WSW_IDU[i] = DrugSexUser['Number of WSW IDU'] / float( DrugUser['Number of IDU agents'])

    HM_NIDU[i] = DrugSexUser['Number of HM NIDU'] / float( DrugUser['Number of NIDU agents'])
    HF_NIDU[i] = DrugSexUser['Number of HF NIDU'] / float( DrugUser['Number of NIDU agents'])
    MSM_NIDU[i] = DrugSexUser['Number of MSM NIDU'] / float( DrugUser['Number of NIDU agents'])
    WSW_NIDU[i] = DrugSexUser['Number of WSW NIDU'] / float( DrugUser['Number of NIDU agents'])

    HM_ND[i] = DrugSexUser['Number of HM ND'] / float( DrugUser['Number of ND agents'])
    HF_ND[i] = DrugSexUser['Number of HF ND'] / float( DrugUser['Number of ND agents'])
    MSM_ND[i] = DrugSexUser['Number of MSM ND'] / float( DrugUser['Number of ND agents'])
    WSW_ND[i] = DrugSexUser['Number of WSW ND'] / float( DrugUser['Number of ND agents'])

    # Number of HIV cases within IDUs	
    IDU_HIV = myModel.get_info_HIV_IDU()
    IDU_HIV_MSM[i] = IDU_HIV['Number of IDU HIV MSM']/float(DrugSexUser['Number of MSM IDU'])
    IDU_HIV_HM[i] = IDU_HIV['Number of IDU HIV HM'] /float(DrugSexUser['Number of HM IDU'])
    IDU_HIV_WSW[i] = IDU_HIV['Number of IDU HIV WSW']/float(DrugSexUser['Number of WSW IDU'])
    IDU_HIV_HF[i] = IDU_HIV['Number of IDU HIV HF']/float(DrugSexUser['Number of HF IDU'])


print('\nDrug User Prevalences:')
print('    IDU          %f +/- %f' % (IDU.mean(), IDU.std()))
print('    NIDU	        %f +/- %f' % (NIDU.mean(), NIDU.std()))
print('    ND           %f +/- %f' % (ND.mean(), ND.std()))

print('\nSex behavior Prevalences:')
print('    HM		%f +/- %f' % (HM.mean(), HM.std()))
print('    HF		%f +/- %f' % (HF.mean(), HF.std()))
print('    MSM		%f +/- %f' % (MSM.mean(), MSM.std()))
print('    WSW		%f +/- %f' % (WSW.mean(), WSW.std()))

print('\nDrug behavior conditioned by Sex behavior:')
print('    IDU  | MSM		%f +/- %f' % (MSM_IDU.mean(), MSM_IDU.std()))
print('    IDU  | HM		%f +/- %f' % (HM_IDU.mean(), HM_IDU.std()))
print('    IDU  | WSW		%f +/- %f' % (WSW_IDU.mean(), WSW_IDU.std()))
print('    IDU  | HF		%f +/- %f' % (HF_IDU.mean(), HF_IDU.std()))

print('    NIDU | MSM		%f +/- %f' % (MSM_NIDU.mean(), MSM_NIDU.std()))
print('    NIDU | HM		%f +/- %f' % (HM_NIDU.mean(), HM_NIDU.std()))
print('    NIDU | WSW		%f +/- %f' % (WSW_NIDU.mean(), WSW_NIDU.std()))
print('    NIDU | HF		%f +/- %f' % (HF_NIDU.mean(), HF_NIDU.std()))

print('    ND   | MSM		%f +/- %f' % (MSM_ND.mean(), MSM_ND.std()))
print('    ND   | HM		%f +/- %f' % (HM_ND.mean(), HM_ND.std()))
print('    ND   | WSW		%f +/- %f' % (WSW_ND.mean(), WSW_ND.std()))
print('    ND   | HF		%f +/- %f' % (HF_ND.mean(), HF_ND.std()))


print('\nPrevalence of HIV within IDUs:')
print('    IDU HIV MSM		%f +/- %f' % (IDU_HIV_MSM.mean(), IDU_HIV_MSM.std()))
print('    IDU HIV HM		%f +/- %f' % (IDU_HIV_HM.mean(), IDU_HIV_HM.std()))
print('    IDU HIV HF		%f +/- %f' % (IDU_HIV_HF.mean(), IDU_HIV_HF.std()))
print('    IDU HIV WSW		%f +/- %f' % (IDU_HIV_WSW.mean(), IDU_HIV_WSW.std()))


