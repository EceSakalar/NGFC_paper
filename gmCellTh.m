function [ output_args ] = gmCellTh(bs, ch_n)
%This program analyses the gamma oscillatory coupling relative to CSD
%traces (spectral aproach; multichannel LFP and CSD traces)
for je=1:ch_n
    for un_indx=1:length(un_chans)
        gm_cellcwtspect{je,un_indx}=[];
    end
    clear un_indx
end
clear je

for je=1:ch_n %CSD channel stepper
    fl_name=strcat(bs.name, bs.num, bs.exp, 'WAW.mat'); %file name
    th_pname=strcat(bs.name, bs.num, bs.exp, 'THP.txt'); %theta periods
    
    if exist(th_pname,'file')==2
        th_pers=dlmread(th_pname);
        clear th_pname;
        
        %load the csd traces
        csd_name=strcat(bs.name, bs.num, bs.exp, 'CSD.mat');
        load(csd_name);
        
        %check if it is edge and if it is defected channel
        if ~isnan(CSD(je,1))&&isempty(ch_def(ch_def==(je)))
            %invert the CSD so that the sink is down and the source is up
            CSDchan=CSD(je,:)*-1;
            clear CSD csd_name;
            
            %reading the spike sequence for the cells, transforming it to sample number instead of time
            for un_indx=1:length(un_chans)
                ch_name=strcat(bs.name, bs.typ, bs.num, bs.exp,'_Ch', num2str(un_chans(un_indx)));
                ch_load=load(num2str(fl_name),num2str(ch_name));
                fl_spks{un_indx}=((ch_load.(ch_name).('times')).*ch_rate);
                clear ch_name ch_load;
            end
            
            %calculate the theta segment wavelet transforms
            for ka=1:size(th_pers,1) %theta period stepper
                sg_start=(th_pers(ka,1)*ch_rate);
                sg_end=(th_pers(ka,2)*ch_rate);
                sg_CSD=CSDchan((round(sg_start)-0.5*ch_rate):(round(sg_end)+0.5*ch_rate));
                %calculating the wavelet transform for the actual channel of the actual segment.
                sg_cwt=cwt(sg_CSD, gm_scsta:gm_scint:gm_scfin, gm_wavelet{1});
                sg_cwt=conj(sg_cwt);
                sg_cwt=sg_cwt(:,(1+0.5*ch_rate):(1+0.5*ch_rate+ceil(sg_end)-floor(sg_start)));
                
                %extract spike triggered CWT for the particular segment for all units and add to a 2D cell array
                for un_indx=1:length(un_chans)
                    sg_spks=round((fl_spks{un_indx}(fl_spks{un_indx}>sg_start&fl_spks{un_indx}<sg_end)-sg_start+1)');
                    sg_gspect=sg_cwt(:,sg_spks);
                    if isempty(gm_cellcwtspect{je,un_indx})
                        gm_cellcwtspect{je,un_indx}=sg_gspect;
                    else
                        gm_cellcwtspect{je,un_indx}=[gm_cellcwtspect{je,un_indx} sg_gspect];
                    end
                    clear sg_spks sg_gspect
                end 
                clear un_indx
                clear sg_cwt sg_start sg_end sg_CSD
            end 
            clear ka th_pers CSDchan
            strcat('CSD_channel_', num2str(je), '_done for all units')
            clear fl_spks
        else
            strcat('CSD_channel_', num2str(je), '_is_not included')
        end
    end
    
end 
clear je
end

