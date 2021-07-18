# Hopf_consciousness



Codes used in the paper titled "Loss of consciousness reduces the stability of brain hubs and the heterogeneity of brain dynamics".

The goal of this paper is to study the fMRI brain dynamics and their mechanism in patients that suffered brain injuries leading to a disorder of consciousness and from healthy subjects undergoing propofol-induced sedation. The paper consists on three parts; 

1) Empirical measures that describe the fMRI brain dynamics described in the phase-locking matrices calculated using the phases of the brain ROIs given by the Hilbert transform. The phase-interaction matrices of the patients can be found in EBRAINS (link). We use measures of integration, segregation, phase-interaction fluctuations and functional connectivity dynamics. The details of each measure can be found in the original paper. Script: phase_synchronization_measures.m


2) Whole-brain modelling:we used a whole-brain model based on Hopf bifurcations. This model combines single-node local oscillatory dynamics and network interactions. It is able to generate different collective dynamics depending on the shape of anatomical connectivity, the global strength of connections and the local state of the network’s nodes.  Importantly, the model allows investigating the interplay between the network structure and the dynamics at the local and global level. As explained in the original paper, three versions of the modelling were used; 1) homogeneous model, where all nodes were set to a=0 and the goal was to describe the global properties of the spatio-temporal dynamics, 2) heterogeneous model, where we extended the model to allow differences in bifurcation parameters aj for different ROIs (the g parameter was the one estimated with the homogeneous model-they can be found in the Table 2 of the paper-) and 3) we studied the relation between local and network dynamics by changing the heteregoneous model to define the effective bifurcation parameters. The codes for these three models are homogeneous_model.m, heterogeneous_model.m and effective_heterogeneous_model.m, respectively.


3) Graph metric applied in the structural connectivity matrices. The python codes of the structural connectivity analysis are available on (GAlib: Graph Analysis library in Python / Numpy, https://github.com/decolab/pyGAlib).


