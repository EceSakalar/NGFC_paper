function [thcell] = thUntCoupl(LFP, ch)
%Theta unit coupling CALCULUS
fl_name=strcat(bs.name, bs.num, bs.exp, 'WAW.mat');   %file name

%THETA ANALYSIS
%reading pyramidal layer LFP

LFPred=downsample(smooth(LFP(ch.pyr,:),ch.dwnsmpl),ch.dwnsmpl); %downsample the pyramidal layer LFP
clear LFP;

%creating theta filter
th_filter=dfilt.dffir(fir1(flt_ncoeff, [2.*(flt_low./ch.ThRate) 2.*(flt_high./ch.ThRate)],'bandpass', gausswin(flt_ncoeff+1)));

%filtering for theta
LFPflt=filtfilt(th_filter.Numerator,1,LFPred')'; %need transpose or not???? (the other version does not)
clear LFPred;

%reading theta periods
th_pname=strcat(bs.name, bs.num, bs.exp, 'THP.txt');


%detecting theta troughs for pyramidal layer LFP
LFPthreshold=0.001*(mean(abs(LFPflt)));
[~,th_troughsLFP]=findpeaks(LFPflt.*-1,'MINPEAKHEIGHT', LFPthreshold, 'MINPEAKDISTANCE', floor(ch_rate/flt_high)-3);
th_troughsLFP=(th_troughsLFP./ch.ThRate);%-(1/ch_rate)+(1/(ch_rate*ch_dwnsmpl));

for un_indx=1:length(un_chans)
    %reading theta periods
    th_pers=dlmread(th_pname);
    
    %limiting theta periods if necessary
    if un_shnk(un_indx)==0
        for je=1:size(th_pers,1)
            if th_pers(je,2)<fl_limS(1, un.nums(un_indx))
                th_pers(je,:)=0;
            elseif th_pers(je,1)<fl_limS(1, un.nums(un_indx))
                th_pers(je,1)=fl_limS(1, un.nums(un_indx));
            end
            if th_pers(je,1)>fl_limE(1, un.nums(un_indx))
                th_pers(je,:)=0;
            elseif th_pers(je,2)>fl_limE(1, un.nums(un_indx))
                th_pers(je,2)=fl_limE(1, un.nums(un_indx));
            end
        end
        clear je
        th_persN=[th_pers(th_pers(:,1)>0,1) th_pers(th_pers(:,2)>0,2)];
        th_pers=th_persN; clear th_persN
    end
    
    %calculating theta time
    fl_th_time(un_indx)=sum(th_pers(:,2)-th_pers(:,1));
    th_time(un_indx)=th_time(un_indx)+fl_th_time(un_indx);
    
    %reading the spike sequence for the cell
    ch_name=strcat(bs.name, bs.typ, bs.num, bs.exp,'_Ch', num2str(un.chans(un_indx)));
    ch_load=load(num2str(fl_name),num2str(ch_name));
    fl_spks=(ch_load.(ch_name).('times')); clear ch_name ch_load;
    
    fl_thspks=periodcut(fl_spks,th_pers,1);
    
    if isempty(fl_thspks)
        th_coupvect_LFPact=[];
    else
        for el=1:length(fl_thspks)
            spk_act=fl_thspks(el);
            cyc_beg=max(th_troughsLFP(th_troughsLFP<=spk_act));
            cyc_end=min(th_troughsLFP(th_troughsLFP>spk_act));
            cyc_LFP=LFPflt((cyc_beg*ch_rate):(cyc_end*ch_rate));
            [~,cyc_cent]=max(cyc_LFP);
            cyc_cent=(cyc_cent/ch_rate)+cyc_beg;
            if  spk_act<=cyc_cent
                th_coupvect_LFPact(el)=180*((spk_act-cyc_beg)/(cyc_cent-cyc_beg));
            else
                th_coupvect_LFPact(el)=180*((spk_act-cyc_cent)/(cyc_end-cyc_cent))+180;
            end
            clear spk_act cyc_beg cyc_cent cyc_end cyc_LFP;
            
        end
        clear el;
    end
    
    %cycle statistics - spike count
    fl_th_cycspkcount=histc(fl_thspks,th_troughsLFP);
    fl_th_cycspkcount(fl_th_cycspkcount==0)=[];
    %cycle statistics - cycle count
    fl_th_cyccount=length(periodcut(th_troughsLFP,th_pers,1))+size(th_pers,1);
    th_cyccount(un_indx)=th_cyccount(un_indx)+fl_th_cyccount;
    clear fl_th_cyccount
    
    %sum coupling vectors accross files
    if isempty(th_coupvect_LFP{un_indx})
        th_coupvect_LFP{un_indx}=th_coupvect_LFPact;
        th_cycspkcount{un_indx}=fl_th_cycspkcount;
    else
        th_coupvect_LFP{un_indx}=[th_coupvect_LFP{un_indx} th_coupvect_LFPact];
        th_cycspkcount{un_indx}=[th_cycspkcount{un_indx}; fl_th_cycspkcount];
    end
    clear th_coupvect_LFPact fl_th_cycspkcount
    
    %bin time calculus
    for el=1:size(th_pers,1)
        perdur=th_pers(el,2)-th_pers(el,1);
        th_acttroughs=periodcut(th_troughsLFP',th_pers(el,:),1);
        
        %calculates cycle boundaries for each cycle in this THperiod
        for em=1:(length(th_acttroughs)+1)
            if em==1
                cyc_end=th_acttroughs(em);
                cyc_beg=max(th_troughsLFP(th_troughsLFP<th_acttroughs(em)));
            elseif em<=length(th_acttroughs)
                cyc_end=th_acttroughs(em);
                cyc_beg=th_acttroughs(em-1);
            else
                cyc_end=min(th_troughsLFP(th_troughsLFP>th_acttroughs(em-1)));
                cyc_beg=th_acttroughs(em-1);
            end
            cyc_LFP=LFPflt((cyc_beg*ch_rate):(cyc_end*ch_rate));
            [~,cyc_cent]=max(cyc_LFP);
            cyc_cent=(cyc_cent/ch_rate)+cyc_beg;
            cycs(em,1)=cyc_beg;
            cycs(em,2)=cyc_cent;
            cycs(em,3)=cyc_end;
            clear cyc_beg cyc_end cyc_cent;
        end
        clear em;
        
        cycs=cycs-th_pers(el,1);
        cycdur_p=zeros(size(cycs,1),2);
        
        %cycle calculates cycle half lengths within a period
        for em=1:size(cycs,1)
            if em==1
                if cycs(em,2)<=0
                    cycdur_p(em,2)=cycs(em,3);
                else
                    cycdur_p(em,1)=cycs(em,2);
                    cycdur_p(em,2)=cycs(em,3)-cycs(em,2);
                end
            elseif em==size(cycs,1)
                if cycs(em,2)>perdur
                    cycdur_p(em,1)=perdur-cycs(em,1);
                else
                    cycdur_p(em,1)=cycs(em,2)-cycs(em,1);
                    cycdur_p(em,2)=perdur-cycs(em,2);
                end
            else
                cycdur_p(em,1)=cycs(em,2)-cycs(em,1);
                cycdur_p(em,2)=cycs(em,3)-cycs(em,2);
            end
        end
        
        %summing the period data
        if isempty (cycdur{un_indx})
            cycdur{un_indx}=cycdur_p;
        else
            cycdur{un_indx}=[cycdur{un_indx}; cycdur_p];
        end
        clear cycdur_p em cycs perdur
    end
    clear el
end
clear un_indx
clear th_pname

%%%%%%%%%%%%

%Calculating coupling for all channels for the cell
%LFP
binscale=(th_binsize/2):(th_binsize):(360-(th_binsize/2));
for un_indx=1:length(un_chans)
    if isempty(th_coupvect_LFP{un_indx})
        th_coupvect_LFPm(un_indx)=NaN;
        th_coupvect_LFPr(un_indx)=NaN;
        th_coupvect_LFPh(un_indx,:)=NaN;
        th_coupvect_LFPp(un_indx)=NaN;
    else
        th_coupvect_LFPm(un_indx)=rad2deg(circ_mean(deg2rad(th_coupvect_LFP{un_indx}')));
        if th_coupvect_LFPm(un_indx)<0
            th_coupvect_LFPm(un_indx)=th_coupvect_LFPm(un_indx)+360;
        end
        th_coupvect_LFPr(un_indx)=circ_r(deg2rad(th_coupvect_LFP{un_indx}'));
        th_coupvect_LFPh(un_indx,:)=hist(rad2deg(deg2rad(th_coupvect_LFP{un_indx}')),binscale);
        n=length(th_coupvect_LFP{un_indx});
        r=th_coupvect_LFPr(un_indx);
        eR=n*r;
        z=(eR.^2)/n;
        P=exp((sqrt(1+(4*n)+(4*(n.^2-eR.^2))))-(1+(2*n)));
        th_coupvect_LFPp(un_indx)=P;
        clear n P r z eR;
    end
end
clear un_indx;

%creating saving variables
thcell.units=[un_nums; un_shnk];

thcell.cyccount=th_cyccount;
thcell.spikespercycles=th_cycspkcount;
thcell.pyrlayer=ch_pyr;
thcell.smplrate=ch_rate;
thcell.cycbintime=cycbintimes;

thcell.coupvect_lfp=th_coupvect_LFP;
thcell.coupvect_lfp_r=th_coupvect_LFPr;
thcell.coupvect_lfp_m=th_coupvect_LFPm;
thcell.coupvect_lfp_hr=th_coupvect_LFPh./repmat(cycbintimes,(size(th_coupvect_LFPh,1)),1);
thcell.coupvect_lfp_h=th_coupvect_LFPh;
thcell.coupvect_lfp_p=th_coupvect_LFPp;
end

