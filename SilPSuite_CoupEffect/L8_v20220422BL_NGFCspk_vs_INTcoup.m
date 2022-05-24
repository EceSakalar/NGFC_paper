%This program reads spike trains for given NGFCs and also individual INTs
%and calculates the INT cell coupling to mid-gamma oscillations before and
%after NGFC firing
%NGFC spike inclusion criterion: previous spike not closer than 51,5 ms
%(gSyn declines below 0.35)

clear all

%EXPERIMENT identification data
bs.name=['B']; bs.num=['193']; bs.exp=['b']; bs.ngfc=['s2u33'];
an.int=[1]; %give the INT units to be analysed in the current session 
intname='glass01'; %give the int name

an.ngfc=[14]; %give the number (sorrend szam!) for the NGFC

fl_name=strcat(bs.name,bs.num,bs.exp,'_gammaunits_CA1mid.mat');
%load the file for the appropriate gamma analysis (only the spks variable)
load(fl_name,'spks'); clear fl_name

%NGFC spikes
spks_ngf=spks{an.ngfc};

%INT spikes
spks_int=spks{an.int};
clear spks

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
gmCSD_ts=gmCSD_pts(gmCSD_pts(:,3)<0,:);
gmCSD_ps=gmCSD_pts(gmCSD_pts(:,3)>0,:);

%cycling across the NGFC spikes
spks_ngfsel=[];
for ka=1:7; spks_intcyc{ka}=[];end; clear ka;
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
        pks_arnd=([pks_bef(end-3:end,:);pks_aft(1:4,:)]);
        clear pks_bef pks_aft
        for ka=1:7
            spks_intcurr=spks_int(spks_int(:,1)>(pks_arnd(ka)./gmCSD_rate)&spks_int(:,1)<=(pks_arnd(ka+1)./gmCSD_rate),:);
            if isempty(spks_intcyc{ka})
                spks_intcyc{ka}=spks_intcurr;
            else
                spks_intcyc{ka}=[spks_intcyc{ka};spks_intcurr];
            end
            clear spks_intcurr
        end; clear ka
        clear pks_arnd cntr
    end 
    clear isi
end; clear je
        
for ka=1:7
    spks_intcyc{2,ka}=spks_intcyc{1,ka}((spks_intcyc{1,ka}(:,8))==1,:);
end; clear ka

for ka=1:7
    intstat.Rsel(ka)=circ_r(deg2rad(spks_intcyc{2,ka}(:,9)));
    intstat.Ssel(ka)=circ_rtest(deg2rad(spks_intcyc{2,ka}(:,9)));
    intstat.Nsel(ka)=length(deg2rad(spks_intcyc{2,ka}(:,9)));
    intstat.Msel(ka)=rad2deg(circ_mean(deg2rad(spks_intcyc{2,ka}(:,9))));
    intstat.Rall(ka)=circ_r(deg2rad(spks_intcyc{1,ka}(:,9)));
    intstat.Sall(ka)=circ_rtest(deg2rad(spks_intcyc{1,ka}(:,9)));
    intstat.Nall(ka)=length(deg2rad(spks_intcyc{1,ka}(:,9)));
    intstat.Mall(ka)=rad2deg(circ_mean(deg2rad(spks_intcyc{1,ka}(:,9))));
end
clear ka gmCSD_ps gmCSD_pts gmCSD_ts

%saving the relevant data
datum = date;
writename=strcat(bs.name, bs.num, bs.exp,bs.ngfc,'_gCA1mid_int',intname,'coup_NGFCtriggered' );
save(writename, 'bs', 'an', 'gmCSD_rate', 'intstat', 'spks_ngf', 'spks_int', 'spks_ngfsel', 'spks_intcyc', 'datum');

masolas=[intstat.Rsel';intstat.Ssel';intstat.Msel';intstat.Nsel'] 

        

