.. _documentation:

********************
Documentation
********************

.. module:: HIVABM_Evolution
   :synopsis: Population including social network.

Inheritance Diagram
===================

.. py:currentmodule:: HIVABM_Evolution
.. inheritance-diagram:: HIVABM_Evolution


Description
=============

This module implements the :class:`SocialNetworkClass` class which
models the population and the associated social network in out HIV model. 
It contains the following classes:

.. py:currentmodule:: HIVABM_Population

+ :class:`PopulationClass` -- Class describing the population.

.. py:currentmodule:: HIVABM_Network

+ :class:`SocialNetworkClass` -- Class modeling the social network.

.. py:currentmodule:: HIVABM_Evolution

+ :class:`HIVModel` -- Class modeling the evolution of our poplulation.


Reference
===================


PopulationClass
----------------------------

.. py:currentmodule:: HIVABM_Population
.. autoclass:: PopulationClass
   :members:
   :show-inheritance:
   :undoc-members:




SocialNetworkClass
------------------------------

.. py:currentmodule:: HIVABM_Network
.. autoclass:: SocialNetworkClass
   :members:
   :show-inheritance:
   :undoc-members:


Evolution
------------------------------

.. module:: HIVABM_Evolution
.. autoclass:: HIVModel
   :members: _update_population, _needle_transmission,
      _sex_transmission, _drug_transition, _update_IDU, _update_NIDU_ND,
      _VCT, _SEP, _drug_treatment,  _initiate_HAART, _update_AllAgents, run, store_results,
      get_HIV_prevalence_drugs, get_HIV_prevalence, plot_results
   :show-inheritance:
   :undoc-members:



