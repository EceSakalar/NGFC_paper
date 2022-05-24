%This program reads spike trains for pooled pyramidal cells
%and calculates the PYR cell coupling to mid-gamma oscillations before and
%after the maximum amplitude of mid-gamma

clear all

load('C:\_work_esakalar\Summary_Analysis\Sakalar et al., 2021\Report\Revision\gmAmpl_NGFC\maxGammaAmplitude_NGFC.mat')
load('C:\_work_esakalar\Summary_Analysis\Sakalar et al., 2021\NGF_spikes.mat')
load('C:\_work_esakalar\Summary_Analysis\Sakalar et al., 2021\NGF_gm_stats.mat')
cellID=readtable('C:\_work_esakalar\Summary_Analysis\NGF_summary\NGF_cells.xlsx'); %for the identitiy of the cells

%EXPERIMENT identification data
bs.name=['B']; bs.num=['182']; bs.exp=['a'];bs.ngfc=['s0u22'];
% an.nums=[2,6,7,8,19,22,84,31,38,41,48,56,57,58];
% an.units=find(ismember(un.nums, an.nums)); %give the PYR units included in the current session 
%an.units=[5,6,7,8,9,10,11,12,13,14,15,16,17,18];
an.chans=[25,27,28,30,31,33,34,35,36,37,38,39,41,43,45,47,48,49,50,51,52,53,55];
an.units=find(ismember(un.chan, an.chans));

fl_name=strcat(bs.name,bs.num,bs.exp,'_gammaunits_CA1mid.mat');
%load the file for the appropriate gamma analysis (only the spks variable)
load(fl_name,'spks'); clear fl_name

%pool pyr cell gamma data
spks_pyr=[];
for je=1:length(an.units)
    actunit=an.units(je);
    if isempty(spks_pyr)
        spks_pyr=spks{actunit};
    else
        spks_pyr=[spks_pyr; spks{actunit}];
    end
end; clear je actunit 
[~,indx]=sort(spks_pyr(:,1));
spks_pyrsrt=spks_pyr(indx,:);
clear indx spks_pyr

%import the gamma peaks&troughs 
fl_name=strcat(bs.name,bs.num,bs.exp,'_gammacomponent_CA1mid.mat');
%load the file for the appropriate gamma analysis (only the spks variable)
load(fl_name,'output');
load(fl_name,'input');
clear fl_name
gmCSD_pts=output.GMcycle';
gmCSD_rate=input.gamma.rate;
clear input output

%extract the gamma cycle troughs
gmCSD_ts=gmCSD_pts(gmCSD_pts(:,3)<0,:); %troughs
gmCSD_ps=gmCSD_pts(gmCSD_pts(:,3)>0,:); %peaks

% trough_times_wNGFC=gmAmpl_NGFC{14}(gmAmpl_NGFC{14}(:,7)==1&gmAmpl_NGFC{15}(:,7)==1,2);
% trough_times_woNGFC=gmAmpl_NGFC{14}(gmAmpl_NGFC{14}(:,7)==0&gmAmpl_NGFC{15}(:,7)==0,2);
% trough_times_NGFC15=gmAmpl_NGFC{14}(gmAmpl_NGFC{14}(:,7)==1&gmAmpl_NGFC{15}(:,7)==0,2);
% trough_times_NGFC16=gmAmpl_NGFC{14}(gmAmpl_NGFC{14}(:,7)==0&gmAmpl_NGFC{15}(:,7)==1,2);

trough_times_wNGFC=gmAmpl_NGFC{12}(gmAmpl_NGFC{12}(:,7)==1,2);
trough_times_woNGFC=gmAmpl_NGFC{12}(gmAmpl_NGFC{12}(:,7)==0,2);

for ka=1:7; spks_pyrcyc{ka}=[];end;clear ka;

%% Pyr gm coupling for max gm with both nGFC spikes

for je=5:size(trough_times_wNGFC,1)-4
 
    %extracting the pyramidal spikes 
    cntr=trough_times_wNGFC(je,1);
    if cntr~=0
        pks_bef=gmCSD_ps((gmCSD_ps(:,1)<=cntr),:);
        pks_aft=gmCSD_ps((gmCSD_ps(:,1)>cntr),:);
        pks_arnd=([pks_bef(end-3:end,:);pks_aft(1:4,:)]);
        clear pks_bef pks_aft

        for ka=1:7
                spks_pyrcurr=spks_pyrsrt(spks_pyrsrt(:,1)>(pks_arnd(ka)./gmCSD_rate)&spks_pyrsrt(:,1)<=(pks_arnd(ka+1)./gmCSD_rate),:);
                if isempty(spks_pyrcyc{ka})
                    spks_pyrcyc{ka}=spks_pyrcurr;
                else
                    spks_pyrcyc{ka}=[spks_pyrcyc{ka};spks_pyrcurr];
                end
                clear spks_pyrcurr
        end; clear ka
        clear pks_arnd cntr
    end
end 

        
for ka=1:7
    spks_pyrcyc{2,ka}=spks_pyrcyc{1,ka}((spks_pyrcyc{1,ka}(:,8))==1,:);
end; clear ka

for ka=1:7
    pyrstat_wNGFC.Rsel(ka)=circ_r(deg2rad(spks_pyrcyc{2,ka}(:,9)));
    pyrstat_wNGFC.Ssel(ka)=circ_rtest(deg2rad(spks_pyrcyc{2,ka}(:,9)));
    pyrstat_wNGFC.Nsel(ka)=length(deg2rad(spks_pyrcyc{2,ka}(:,9)));
    pyrstat_wNGFC.Msel(ka)=rad2deg(circ_mean(deg2rad(spks_pyrcyc{2,ka}(:,9))));
    pyrstat_wNGFC.Rall(ka)=circ_r(deg2rad(spks_pyrcyc{1,ka}(:,9)));
    pyrstat_wNGFC.Sall(ka)=circ_rtest(deg2rad(spks_pyrcyc{1,ka}(:,9)));
    pyrstat_wNGFC.Nall(ka)=length(deg2rad(spks_pyrcyc{1,ka}(:,9)));
    pyrstat_wNGFC.Mall(ka)=rad2deg(circ_mean(deg2rad(spks_pyrcyc{1,ka}(:,9))));
end
clear ka spks_pyrcyc
for ka=1:7; spks_pyrcyc{ka}=[];end;clear ka;


%% Pyr gm coupling for max gm without any nGFC spikes

for je=5:size(trough_times_woNGFC,1)-4
 
    %extracting the pyramidal spikes 
    cntr=trough_times_woNGFC(je,1);
    if cntr~=0
        pks_bef=gmCSD_ps((gmCSD_ps(:,1)<=cntr),:);
        pks_aft=gmCSD_ps((gmCSD_ps(:,1)>cntr),:);
        pks_arnd=([pks_bef(end-3:end,:);pks_aft(1:4,:)]);
        clear pks_bef pks_aft

        for ka=1:7
                spks_pyrcurr=spks_pyrsrt(spks_pyrsrt(:,1)>(pks_arnd(ka)./gmCSD_rate)&spks_pyrsrt(:,1)<=(pks_arnd(ka+1)./gmCSD_rate),:);
                if isempty(spks_pyrcyc{ka})
                    spks_pyrcyc{ka}=spks_pyrcurr;
                else
                    spks_pyrcyc{ka}=[spks_pyrcyc{ka};spks_pyrcurr];
                end
                clear spks_pyrcurr
        end; clear ka
        clear pks_arnd cntr
    end
end 

        
for ka=1:7
    spks_pyrcyc{2,ka}=spks_pyrcyc{1,ka}((spks_pyrcyc{1,ka}(:,8))==1,:);
end; clear ka

for ka=1:7
    pyrstat_woNGFC.Rsel(ka)=circ_r(deg2rad(spks_pyrcyc{2,ka}(:,9)));
    pyrstat_woNGFC.Ssel(ka)=circ_rtest(deg2rad(spks_pyrcyc{2,ka}(:,9)));
    pyrstat_woNGFC.Nsel(ka)=length(deg2rad(spks_pyrcyc{2,ka}(:,9)));
    pyrstat_woNGFC.Msel(ka)=rad2deg(circ_mean(deg2rad(spks_pyrcyc{2,ka}(:,9))));
    pyrstat_woNGFC.Rall(ka)=circ_r(deg2rad(spks_pyrcyc{1,ka}(:,9)));
    pyrstat_woNGFC.Sall(ka)=circ_rtest(deg2rad(spks_pyrcyc{1,ka}(:,9)));
    pyrstat_woNGFC.Nall(ka)=length(deg2rad(spks_pyrcyc{1,ka}(:,9)));
    pyrstat_woNGFC.Mall(ka)=rad2deg(circ_mean(deg2rad(spks_pyrcyc{1,ka}(:,9))));
end
clear ka spks_pyrcyc
for ka=1:7; spks_pyrcyc{ka}=[];end;clear ka;



%% Pyr gm coupling for max gm with only nGFC u15 spikes

for je=5:size(trough_times_NGFC15,1)-4
 
    %extracting the pyramidal spikes 
    cntr=trough_times_NGFC15(je,1);
    if cntr~=0
        pks_bef=gmCSD_ps((gmCSD_ps(:,1)<=cntr),:);
        pks_aft=gmCSD_ps((gmCSD_ps(:,1)>cntr),:);
        pks_arnd=([pks_bef(end-3:end,:);pks_aft(1:4,:)]);
        clear pks_bef pks_aft

        for ka=1:7
                spks_pyrcurr=spks_pyrsrt(spks_pyrsrt(:,1)>(pks_arnd(ka)./gmCSD_rate)&spks_pyrsrt(:,1)<=(pks_arnd(ka+1)./gmCSD_rate),:);
                if isempty(spks_pyrcyc{ka})
                    spks_pyrcyc{ka}=spks_pyrcurr;
                else
                    spks_pyrcyc{ka}=[spks_pyrcyc{ka};spks_pyrcurr];
                end
                clear spks_pyrcurr
        end; clear ka
        clear pks_arnd cntr
    end
end 

        
for ka=1:7
    spks_pyrcyc{2,ka}=spks_pyrcyc{1,ka}((spks_pyrcyc{1,ka}(:,8))==1,:);
end; clear ka

for ka=1:7
    pyrstat_NGFC15.Rsel(ka)=circ_r(deg2rad(spks_pyrcyc{2,ka}(:,9)));
    pyrstat_NGFC15.Ssel(ka)=circ_rtest(deg2rad(spks_pyrcyc{2,ka}(:,9)));
    pyrstat_NGFC15.Nsel(ka)=length(deg2rad(spks_pyrcyc{2,ka}(:,9)));
    pyrstat_NGFC15.Msel(ka)=rad2deg(circ_mean(deg2rad(spks_pyrcyc{2,ka}(:,9))));
    pyrstat_NGFC15.Rall(ka)=circ_r(deg2rad(spks_pyrcyc{1,ka}(:,9)));
    pyrstat_NGFC15.Sall(ka)=circ_rtest(deg2rad(spks_pyrcyc{1,ka}(:,9)));
    pyrstat_NGFC15.Nall(ka)=length(deg2rad(spks_pyrcyc{1,ka}(:,9)));
    pyrstat_NGFC15.Mall(ka)=rad2deg(circ_mean(deg2rad(spks_pyrcyc{1,ka}(:,9))));
end
clear ka spks_pyrcyc
for ka=1:7; spks_pyrcyc{ka}=[];end;clear ka; 

%% Pyr gm coupling for max gm with only nGFC u16 spikes

for je=5:size(trough_times_NGFC16,1)-4
 
    %extracting the pyramidal spikes 
    cntr=trough_times_NGFC16(je,1);
    if cntr~=0
        pks_bef=gmCSD_ps((gmCSD_ps(:,1)<=cntr),:);
        pks_aft=gmCSD_ps((gmCSD_ps(:,1)>cntr),:);
        pks_arnd=([pks_bef(end-3:end,:);pks_aft(1:4,:)]);
        clear pks_bef pks_aft

        for ka=1:7
                spks_pyrcurr=spks_pyrsrt(spks_pyrsrt(:,1)>(pks_arnd(ka)./gmCSD_rate)&spks_pyrsrt(:,1)<=(pks_arnd(ka+1)./gmCSD_rate),:);
                if isempty(spks_pyrcyc{ka})
                    spks_pyrcyc{ka}=spks_pyrcurr;
                else
                    spks_pyrcyc{ka}=[spks_pyrcyc{ka};spks_pyrcurr];
                end
                clear spks_pyrcurr
        end; clear ka
        clear pks_arnd cntr
    end
end 

        
for ka=1:7
    spks_pyrcyc{2,ka}=spks_pyrcyc{1,ka}((spks_pyrcyc{1,ka}(:,8))==1,:);
end; clear ka

for ka=1:7
    pyrstat_NGFC16.Rsel(ka)=circ_r(deg2rad(spks_pyrcyc{2,ka}(:,9)));
    pyrstat_NGFC16.Ssel(ka)=circ_rtest(deg2rad(spks_pyrcyc{2,ka}(:,9)));
    pyrstat_NGFC16.Nsel(ka)=length(deg2rad(spks_pyrcyc{2,ka}(:,9)));
    pyrstat_NGFC16.Msel(ka)=rad2deg(circ_mean(deg2rad(spks_pyrcyc{2,ka}(:,9))));
    pyrstat_NGFC16.Rall(ka)=circ_r(deg2rad(spks_pyrcyc{1,ka}(:,9)));
    pyrstat_NGFC16.Sall(ka)=circ_rtest(deg2rad(spks_pyrcyc{1,ka}(:,9)));
    pyrstat_NGFC16.Nall(ka)=length(deg2rad(spks_pyrcyc{1,ka}(:,9)));
    pyrstat_NGFC16.Mall(ka)=rad2deg(circ_mean(deg2rad(spks_pyrcyc{1,ka}(:,9))));
end
clear ka spks_pyrcyc

%saving the relevant data

datum = date;
writename=strcat(bs.name, bs.num, bs.exp,bs.ngfc,'_gCA1mid_PYRcoup_NGFCtriggered' );
save(writename, 'bs', 'an', 'gmCSD_rate', 'pyrstat', 'spks_ngf', 'spks_pyrsrt', 'spks_ngfsel', 'spks_pyrcyc', 'datum');



        

