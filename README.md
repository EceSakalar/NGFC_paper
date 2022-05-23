# NGFC_paper
Code for reproducibility of our paper Sakalar et al., 2022

The scripts in this project includes:
 -Current source density (CSD) analysis calculation from local field potentials (LFP)
 -Sharp wave-associated ripple detection in CSD and LFP
 -Cell firing versus sharp wave-associated ripples
 -Theta oscillation detection through band-pass filtering
 -Cell firing versus theta oscillation
 -Gamma oscillation detection through wavelet transformation
 -Cell firing versus gamma oscillation
 -Gamma oscillation detection through band-pass filtering (for the gamma frequency range detected manually)
 -Cell firing versus gamma oscillation
 


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
