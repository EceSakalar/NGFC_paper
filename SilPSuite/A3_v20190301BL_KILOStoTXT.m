%This program reads Kilosort Phy files and exports spikes of each unit into a single .txt file
close all
clear all

%USER INPUT
%give the experimenter identifier, the experiment number and the experimental day
bs_name='B'; bs_num='216'; bs_day='_b';

%give the unit numbers you want to include
phy_units = [0 4 10 11 12 20 23 47 218, 35 39 40 48 49 236, 54 62 190, 64 65 66, 82 92 103 227, 104 106 115, 125 126 129 132 183 215, 141 144 147 149 150 152 153 154]; 
phy_shank = [repmat([1],1,9),repmat([2],1,6),repmat([3],1,3),repmat([4],1,3),repmat([5],1,4),repmat([6],1,3),repmat([7],1,6),repmat([8],1,8)];
%give the sampling rate in Hz
smpl_rate = 20000;
smpl_int = 1/smpl_rate;

%step along all files
un_clstnums=readNPY('spike_clusters.npy');
un_smpls=readNPY('spike_times.npy');

for ka=1:size(phy_units,2)
    un_spktimes=(double(un_smpls(un_clstnums==phy_units(ka)))).*smpl_int;
    tx_name=strcat(bs_name, bs_num, bs_day, '_s', num2str(phy_shank(ka)),'_u',num2str(phy_units(ka)),'.txt');
    dlmwrite (tx_name,un_spktimes, 'precision', 12);
    clear tx_name un_spktimes
end; clear ka

