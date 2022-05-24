%This program reads spike trains for pooled pyramidal cells
%and calculates the PYR cell coupling to mid-gamma oscillations before and
%after the maximum amplitude of mid-gamma

%clear all

%EXPERIMENT identification data
bs.name=['B']; bs.num=['207']; bs.exp=['b'];bs.ngfc=['s0u22'];
% an.chans=[134	135	136	139	144	145	148	150	153	154	156	157	161	163	164]; 
% an.units=find(ismember(un.chan, an.chans));
% an.shnk=5; %reference shank where ngf would be
% an.units=[1	4	5	6	7	12	13	14	15	22	23	24	26	27	29	31	32	33	34	35];
% fl_name=strcat(bs.name,bs.num,bs.exp,'_gammaunits_CA1mid.mat');
% %load the file for the appropriate gamma analysis (only the spks variable)
% load(fl_name,'spks'); clear fl_name

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
% fl_name=strcat(bs.name,bs.num,bs.exp,'_gammacomponent_CA1mid.mat');
% %load the file for the appropriate gamma analysis (only the spks variable)
% load(fl_name,'output');
% load(fl_name,'input');
% clear fl_name
gmCSD_pts=(output.GMcycle{an.shnk})';
gmCSD_rate=input.gamma.rate;
clear input output

%extract the gamma cycle troughs
gmCSD_ts=gmCSD_pts(gmCSD_pts(:,3)<0,:); %troughs
gmCSD_ps=gmCSD_pts(gmCSD_pts(:,3)>0,:); %peaks

trough_times=gmAmpl_NGFC{25}(:,2);
for ka=1:7; spks_pyrcyc{ka}=[];end;clear ka;

for je=5:size(trough_times,1)-4
 
    %extracting the pyramidal spikes 
    cntr=trough_times(je,1);
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
    pyrstat.Rsel(ka)=circ_r(deg2rad(spks_pyrcyc{2,ka}(:,9)));
    pyrstat.Ssel(ka)=circ_rtest(deg2rad(spks_pyrcyc{2,ka}(:,9)));
    pyrstat.Nsel(ka)=length(deg2rad(spks_pyrcyc{2,ka}(:,9)));
    pyrstat.Msel(ka)=rad2deg(circ_mean(deg2rad(spks_pyrcyc{2,ka}(:,9))));
    pyrstat.Rall(ka)=circ_r(deg2rad(spks_pyrcyc{1,ka}(:,9)));
    pyrstat.Sall(ka)=circ_rtest(deg2rad(spks_pyrcyc{1,ka}(:,9)));
    pyrstat.Nall(ka)=length(deg2rad(spks_pyrcyc{1,ka}(:,9)));
    pyrstat.Mall(ka)=rad2deg(circ_mean(deg2rad(spks_pyrcyc{1,ka}(:,9))));
end
clear ka gmCSD_ps gmCSD_pts gmCSD_ts

%saving the relevant data
% datum = date;
% writename=strcat(bs.name, bs.num, bs.exp,bs.ngfc,'_gCA1mid_PYRcoup_NGFCtriggered' );
% save(writename, 'bs', 'an', 'gmCSD_rate', 'pyrstat', 'spks_ngf', 'spks_pyrsrt', 'spks_ngfsel', 'spks_pyrcyc', 'datum');
% 
% 

        

