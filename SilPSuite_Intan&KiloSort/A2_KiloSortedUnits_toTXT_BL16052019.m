
%USER INPUT
%give the experimenter identifier, the experiment number, the experimental day and the shank number
bs_name='IV'; bs_num='004a'; bs_day='_d1'; 

%give the sampling rate in Hz
smpl_rate=20000;

%read spike times, spike clusters and 
sp_times=readNPY('spike_times.npy');
sp_ident=readNPY('spike_clusters.npy');
un_categ=readtable('cluster_groups.csv');

for je=1:size(un_categ,1)
    if un_categ{je,2}{1}(1)=='g'
        un_actual=un_categ{je,1};
        sp_actual=single(sp_times(sp_ident==un_actual)).*(1/smpl_rate);
        tx_name=strcat(bs_name, bs_num, bs_day, '_u',num2str(un_actual),'.txt');
        dlmwrite (tx_name,sp_actual, 'precision', 12);
        clear tx_name sp_actual un_actual
    end
end