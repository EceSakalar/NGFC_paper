%This program reads .kwik file cl usters and exports spikes of each unit into a single .txt file
close all
clear all

%USER INPUT
%give the experimenter identifier, the experiment number, the experimental day and the shank number
bs_name='B'; bs_num='193'; bs_day='_d2'; bs_shank='_shR';

%give the file names for the .kwik files (comma delimited)
kw_name = {'B193b_R.kwik'};
%file numbers
kw_files = (0);
%give the unit numbers you want to include
kw_units = [20 27 32 70 85 91 98 140 142 147 148 149 158 161 163 164 165 168 169 170]; 

%give the sampling rate in Hz
smpl_rate = 20000;
smpl_int = 1/smpl_rate;

%step along all files
for je=1:size(kw_name,2)
    kw_read=char(kw_name(1,je));
    un_clstnums=hdf5read(kw_read, '/channel_groups/0/spikes/clusters/main');
    un_smpls=hdf5read(kw_read, '/channel_groups/0/spikes/time_samples');
    
    for ka=1:size(kw_units,2)
        un_spktimes=(double(un_smpls(un_clstnums==kw_units(ka)))).*smpl_int;
         
        tx_name=strcat(bs_name, bs_num, bs_day, bs_shank,'_u',num2str(kw_units(ka)),'_f', num2str(kw_files(je)),'.txt');
        dlmwrite (tx_name,un_spktimes, 'precision', 12);
        clear tx_name un_spktimes 
    end; clear ka
    clear kw_read un_smpls un_clstnums
end; clear je