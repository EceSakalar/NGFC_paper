# NGFC_paper
Code for reproducibility of our paper Sakalar et al., 2022

The scripts in this project include:

 -Current source density (CSD) analysis calculation from local field potentials (LFP)
 -Sharp wave-associated ripple detection in CSD and LFP
 -Cell firing versus sharp wave-associated ripples
 -Theta oscillation detection through band-pass filtering
 -Cell firing versus theta oscillation
 -Gamma oscillation detection through wavelet transformation
 -Cell firing versus gamma oscillation
 -Gamma oscillation detection through band-pass filtering (for the gamma frequency range detected manually)
 -Cell firing versus gamma oscillation
 -Gamma oscillation coupling of cells depending on timing of neurogliaform cell firing or gamma amplitude
 -Cell firing versus brain states
 -Synaptic conductance modelling of cell firing (specified on neurogliaform cells)
 -Analysis of position of the animals in the virtual environment
 -Spike waveform characteristics of cells
 -Inter-cell interactions (analysis of inhibition, excitation and co-firing)


Required input:
.mat file with LFP and spike times
.txt files for manually detected theta oscillation and sharp wave-associated ripple times for the LFP in the .mat file

Required system:
MATLAB R2017a
  -Statistics toolbox

References:
For the circular statistics, the functions from CircStat Toolbox was used:
P. Berens, CircStat: A Matlab Toolbox for Circular Statistics, Journal of Statistical Software, Volume 31, Issue 10, 2009
http://www.jstatsoft.org/v31/i10
