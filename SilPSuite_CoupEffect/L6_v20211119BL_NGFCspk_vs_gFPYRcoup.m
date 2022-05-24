%This program reads spike trains for given NGFCs and also pooled pyramidal cells
%and calculates the PYR cell coupling to fast-gamma oscillations before and
%after NGFC firing
%NGFC spike inclusion criterion: previous spike not closer than 51,5 ms
%(gSyn declines below 0.35)

clear all,

%EXPERIMENT identification data
bs.name=['B']; bs.num=['179']; bs.exp=['a'];bs.ngfc=['gls01'];
an.units=[2,4,5,6,7,8,11,12,14,15,16,18,19,20]; %give the PYR units included in the current session 

an.ngfc=[23]; %give the number (sorrend szam!) for the NGFC

fl_name=strcat(bs.name,bs.num,bs.exp,'_gammaunits_CA1fst.mat');
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
%NGFC spikes
spks_ngf=spks{an.ngfc};
clear spks

%import the gamma peaks&troughs 
fl_name=strcat(bs.name,bs.num,bs.exp,'_gammacomponent_CA1fst.mat');
%load the file for the appropriate gamma analysis (only the spks variable)
load(fl_name,'output');
load(fl_name,'input');
clear fl_name
gmCSD_pts=output.GMcycle';
gmCSD_rate=input.gamma.rate;
clear input output

%extract the gamma cycle troughs
gmCSD_ts=gmCSD_pts(gmCSD_pts(:,3)<0,:);
gmCSD_ps=gmCSD_pts(gmCSD_pts(:,3)>0,:);
%cycling across the NGFC spikes
%for gammaM it was 7 cycles (+/-3)
%for gammaF it is 11 cycles (+/-5)
spks_ngfsel=[];
for ka=1:11; spks_pyrcyc{ka}=[];end;clear ka;
for je=1:size(spks_ngf,1)
    if je==1; isi=spks_ngf(1,1); else isi=spks_ngf(je,1)-spks_ngf(je-1,1); end
    %include NGFC spike if isi is >0.0515 s (51.5 ms, where gSyn declines <0.35 after a spike)
    if isi>=0.0515&&spks_ngf(je,3)==1
        if isempty(spks_ngfsel)
            spks_ngfsel=spks_ngf(je,:);
        else
            spks_ngfsel=[spks_ngfsel;spks_ngf(je,:)];
        end
        %extracting the pyramidal spikes 
        cntr=spks_ngf(je,1);
        pks_bef=gmCSD_ps((gmCSD_ps(:,1)<=(cntr.*gmCSD_rate)),:);
        pks_aft=gmCSD_ps((gmCSD_ps(:,1)>(cntr.*gmCSD_rate)),:);
        pks_arnd=([pks_bef(end-5:end,:);pks_aft(1:6,:)]);
        clear pks_bef pks_aft
        for ka=1:11
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
    clear isi
end; clear je
        
for ka=1:11
    spks_pyrcyc{2,ka}=spks_pyrcyc{1,ka}((spks_pyrcyc{1,ka}(:,8))==1,:);
end; clear ka

for ka=1:11
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
datum = date;
writename=strcat(bs.name, bs.num, bs.exp,bs.ngfc,'_gCA1fst_PYRcoup_NGFCtriggered' );
save(writename, 'bs', 'an', 'gmCSD_rate', 'pyrstat', 'spks_ngf', 'spks_pyrsrt', 'spks_ngfsel', 'spks_pyrcyc', 'datum');

sumer=[pyrstat.Rsel';NaN;pyrstat.Ssel';NaN;pyrstat.Msel';NaN;pyrstat.Nsel'];


        

