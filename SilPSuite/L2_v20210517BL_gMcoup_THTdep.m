%This script generates statistics for the Theta modulation of 
%mid-gamma oscillation -COUPLING for single cells or small populations of cells 

clear all

%EXPERIMENT identification data
bs.name=['ES']; bs.num=['33']; bs.exp=['_d3'];
an.units=[2,3,6,10,20,22,23,24]; %give the units included in the current session 

fl_name=strcat(bs.name,bs.num,bs.exp,'_gammaunits_CA1mid.mat');
%load the file for the appropriate gamma analysis (only the spks variable)
load(fl_name,'spks'); clear fl_name


%extract the theta and gamma phases for the included cycles
for je=1:length(an.units)
    actunit=an.units(je);
    spks_actunit=spks{actunit};
    
    el=1;
    for ka=1:size(spks_actunit,1)
        if ~isnan(spks_actunit(ka,3))
            phs_thgm_t(el,1)=spks_actunit(ka,4);
            phs_thgm_t(el,2)=spks_actunit(ka,9);
            el=el+1;
        end
    end; clear ka el
    if je==1
        phs_thgm=phs_thgm_t;
    else
        phs_thgm=[phs_thgm; phs_thgm_t];
    end
    clear actunit phs_thgm_t spks_actunit   
end; clear je

%generate the thta phase binned
%20 bins centered on 9:18:351 (°) each 18° wide (not overlapping bins)
phs_thgmex=[[phs_thgm(:,1)-360,phs_thgm(:,2)];phs_thgm;[phs_thgm(:,1)+360,phs_thgm(:,2)]];
for je=1:20
    lim_L=(je*18)-27; lim_U=(je*18)+9;
    bin_samps{je}=phs_thgmex((phs_thgmex(:,1)<=lim_U&phs_thgmex(:,1)>=lim_L),:);
    
    bin_gMr(je)=circ_r(deg2rad(bin_samps{je}(:,2)));
    bin_gMm(je)=rad2deg(circ_mean(deg2rad(bin_samps{je}(:,2))));
    bin_gMp(je)=circ_rtest(deg2rad(bin_samps{je}(:,2)));
    bin_gMn(je)=size(bin_samps{je},1);
    bin_THTax(je)=lim_L+18;
    clear lim_L lim_U
end; clear je
tocopy=[bin_THTax',bin_gMn',bin_gMr',bin_gMm',bin_gMp'];
for je=1:20
    if tocopy(je,5)<=0.05
       tocopy(je,6)=tocopy(je,3);
       if tocopy(je,4)<=-90
           tocopy(je,7)=tocopy(je,4)+360;
       else
           tocopy(je,7)=tocopy(je,4);
       end
    else
        tocopy(je,6)=NaN;
        tocopy(je,7)=NaN;
    end
end
        

