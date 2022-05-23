%This program generates theta modulation of unit firing relative to 'pyramidal layer' LFP 
clear all
%close all

bs_name='B'; bs_num='183'; bs_exp='a'; bs_typ='sx'; %basic data entry
ch_dwnsmpl=5;                                       %give the downsampling rate
ch_rate=2000; ch_rate=ch_rate/ch_dwnsmpl;           %give the silicon probe sampling rate (in Hz)

un_chans=[23:50,22];                                                           %channels containing the spike times for the units
un_nums=[0,4,5,8,13,20,24,26,27,31,62,72,80,38,44,46,49,55,69,103,109,207,227,228,239,264,280,294,1];          %unit numbers (in kwik file)
%shanks from which unit was isolated; shank 0 is glass electrode and allows shortening of the recording time
un_shnk =[repmat([1],1,13) repmat([2],1,0) repmat([3],1,15) 0];                              

ch_pyr=14;                                       %channel number (ventralmost LFP contact is 1) for the pyramidal layer
fl_ord=[1];                                     %file extension numbers belonging to the present cell

%here you can give the limits of inclusion for each file (rows) for each unit (colums) on glass electrode
fl_limS = [2682]; %starts
fl_limE = [2931]; %ends

%number of theta bins
th_nbins = 20;                                      
th_binsize = 360/th_nbins;

%filter properties 
flt_ncoeff=512; flt_low=5; flt_high=12;

%generating empty variables
for un_indx=1:length(un_chans)
    th_coupvect_LFP{un_indx}=[];
    cycdur{un_indx}=[];
end
th_time=zeros(length(un_nums),1)';
th_cyccount=zeros(length(un_nums),1)';

%CALCULUS
for fl_indx=1:length(fl_ord); %FILE STEPPER
    fl_name=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'WAW.mat');   %file name
    
    %THETA ANALYSIS
    %reading pyramidal layer LFP
    lfp_name=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'LFP.mat');
    load(lfp_name);
    LFPred=downsample(smooth(LFP(ch_pyr,:),ch_dwnsmpl),ch_dwnsmpl);
    clear LFP;
    
    %creating theta filter
    th_filter=dfilt.dffir(fir1(flt_ncoeff, [2.*(flt_low./ch_rate) 2.*(flt_high./ch_rate)],'bandpass', gausswin(flt_ncoeff+1)));
    
    %filtering for theta
    LFPflt=filtfilt(th_filter.Numerator,1,LFPred')'; clear LFPred;
    
    %reading theta periods
    th_pname=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'THP.txt');
    
    if exist(th_pname,'file')==2;
        
        %detecting theta troughs for pyramidal layer LFP
        LFPthreshold=0.001*(mean(abs(LFPflt)));
        [~,th_troughsLFP]=findpeaks(LFPflt.*-1,'MINPEAKHEIGHT', LFPthreshold, 'MINPEAKDISTANCE', floor(ch_rate/flt_high)-3);
        th_troughsLFP=(th_troughsLFP./ch_rate);%-(1/ch_rate)+(1/(ch_rate*ch_dwnsmpl));
        
        for un_indx=1:length(un_chans)
            %reading theta periods
            th_pers=dlmread(th_pname);
            
            %limiting theta periods if necessary
            if un_shnk(un_indx)==0;
                for je=1:size(th_pers,1)
                    if th_pers(je,2)<fl_limS(fl_indx, un_nums(un_indx))
                        th_pers(je,:)=0;
                    elseif th_pers(je,1)<fl_limS(fl_indx, un_nums(un_indx))
                        th_pers(je,1)=fl_limS(fl_indx, un_nums(un_indx));
                    end
                    if th_pers(je,1)>fl_limE(fl_indx, un_nums(un_indx))
                        th_pers(je,:)=0;
                    elseif th_pers(je,2)>fl_limE(fl_indx, un_nums(un_indx))
                        th_pers(je,2)=fl_limE(fl_indx, un_nums(un_indx));
                    end
                end; clear je
                th_persN=[th_pers(th_pers(:,1)>0,1) th_pers(th_pers(:,2)>0,2)];
                th_pers=th_persN; clear th_persN
            end
            
            %calculating theta time
            fl_th_time(un_indx)=sum(th_pers(:,2)-th_pers(:,1));
            th_time(un_indx)=th_time(un_indx)+fl_th_time(un_indx);
            
            %reading the spike sequence for the cell
            ch_name=strcat(bs_name, bs_typ, bs_num, bs_exp,'_Ch', num2str(un_chans(un_indx)));
            ch_load=load(num2str(fl_name),num2str(ch_name));
            fl_spks=(ch_load.(ch_name).('times')); clear ch_name ch_load;
            
            fl_thspks=periodcut(fl_spks,th_pers,1);
            
            if isempty(fl_thspks)
                th_coupvect_LFPact=[];
            else
                for el=1:length(fl_thspks)
                    spk_act=fl_thspks(el);
                    cyc_beg=max(th_troughsLFP(th_troughsLFP<=spk_act));cyc_end=min(th_troughsLFP(th_troughsLFP>spk_act));
                    cyc_LFP=LFPflt((cyc_beg*ch_rate):(cyc_end*ch_rate)); [~,cyc_cent]=max(cyc_LFP);cyc_cent=(cyc_cent/ch_rate)+cyc_beg;
                    if  spk_act<=cyc_cent;      th_coupvect_LFPact(el)=180*((spk_act-cyc_beg)/(cyc_cent-cyc_beg));
                    else                        th_coupvect_LFPact(el)=180*((spk_act-cyc_cent)/(cyc_end-cyc_cent))+180;
                    end; clear spk_act cyc_beg cyc_cent cyc_end cyc_LFP;
                    
                end; clear el;
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
                th_coupvect_LFP{un_indx}=[th_coupvect_LFP{un_indx}  th_coupvect_LFPact];
                th_cycspkcount{un_indx}=[th_cycspkcount{un_indx}; fl_th_cycspkcount];
            end; clear th_coupvect_LFPact fl_th_cycspkcount
            
            %bin time calculus
            for el=1:size(th_pers,1)
                perdur=th_pers(el,2)-th_pers(el,1);
                th_acttroughs=periodcut(th_troughsLFP',th_pers(el,:),1);
                
                %calculates cycle boundaries for each cycle in this THperiod
                for em=1:(length(th_acttroughs)+1);
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
                    [~,cyc_cent]=max(cyc_LFP);cyc_cent=(cyc_cent/ch_rate)+cyc_beg;
                    cycs(em,1)=cyc_beg; cycs(em,2)=cyc_cent;cycs(em,3)=cyc_end;
                    clear cyc_beg cyc_end cyc_cent;
                end; clear em;
                
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
                end;
                
                %summing the period data
                if isempty (cycdur{un_indx})
                    cycdur{un_indx}=cycdur_p;
                else cycdur{un_indx}=[cycdur{un_indx}; cycdur_p];
                end; clear cycdur_p em cycs perdur
            end; clear el
        end; clear un_indx
    end;
end; clear th_pname

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
        end;
        th_coupvect_LFPr(un_indx)=circ_r(deg2rad(th_coupvect_LFP{un_indx}'));
        th_coupvect_LFPh(un_indx,:)=hist(rad2deg(deg2rad(th_coupvect_LFP{un_indx}')),binscale);
        n=length(th_coupvect_LFP{un_indx});
        r=th_coupvect_LFPr(un_indx);
        eR=n*r;
        z=(eR.^2)/n; P=exp((sqrt(1+(4*n)+(4*(n.^2-eR.^2))))-(1+(2*n)));
        masolashoz(1:5,un_indx)=[n; n/th_time(un_indx); th_coupvect_LFPm(un_indx); r; P];
        th_coupvect_LFPp(un_indx)=P; clear n P r z eR;
    end
end; clear un_indx;


%PLOTTING STATISTICS
%plotting figures
figindex=0;

for un_indx=1:length(un_chans)
    figindex_T=ceil(un_indx/6);
    if figindex_T~=figindex
        figindex=figindex+1;
        figname=strcat(bs_name, bs_num, bs_exp, num2str(figindex),'_thetamodulation');
        figure ('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
        plotindex=1;
    else plotindex=((un_indx-((figindex-1)*6))*3)-2;
    end
    if isempty(th_coupvect_LFP{un_indx})
        subplot(6,3,plotindex);
        axis off;
        text (0, 1, 'no data for this cell');
    else
        %textual results
        subplot (6,3,plotindex);
        axis off;
        actualtext= strcat(bs_name, bs_num, bs_exp, 'sh', num2str(un_shnk(un_indx)),'u',num2str(un_nums(un_indx))); text (0, 1, actualtext);
        %number of theta cycles
        actualtext= strcat('Number of theta cycles analysed: ', num2str (th_cyccount(un_indx))); text (0, 0.9, actualtext);
        %number of active theta cycles
        th_acyccount=length(th_cycspkcount{un_indx});th_acycperc=(th_acyccount/th_cyccount(un_indx))*100; actualtext= strcat('of these active: ', num2str (th_acyccount), ' (', num2str(th_acycperc), '%)'); text (0, 0.8, actualtext);
        %%mean number of APs in active theta cycles
        actualtext= strcat('# of spikes/cycle (SD):', num2str(mean(th_cycspkcount{un_indx})), ' (', num2str(std(th_cycspkcount{un_indx})), ')'); text (0, 0.7, actualtext);
        %mean AP frequency
        th_spkrate=length(th_coupvect_LFP{un_indx})/th_time(un_indx); actualtext=strcat('spike frequency under theta:', num2str (th_spkrate), ' s-1'); text (0, 0.6, actualtext);
        %mean vector legth
        actualtext= strcat('mean vector length:', num2str (th_coupvect_LFPr(un_indx))); text (0, 0.5, actualtext);
        %mean phase angle
        actualtext= strcat('mean vector angle:', num2str (th_coupvect_LFPm(un_indx))); text (0, 0.4, actualtext);
        % angular deviation and circular SD
        s=(180/pi)*(sqrt(2*(1-th_coupvect_LFPr(un_indx)))); so=(180/pi)*(sqrt(-2*(log(th_coupvect_LFPr(un_indx)))));
        actualtext=strcat('Angular deviation (s): ', num2str(s)); text (0, 0.3, actualtext);
        actualtext=strcat('Circular standard deviation (so): ', num2str(so)); text (0, 0.2, actualtext);
        %Rayleigh statistics
        n=length(th_coupvect_LFP{un_indx});
        R=n*th_coupvect_LFPr(un_indx); z=R.^2/n; P=exp((sqrt(1+(4*n)+(4*(n.^2-R.^2))))-(1+(2*n)));
        text (0, 0, 'Rayleigh test for the uniformity of distribution');actualtext=strcat ('z=', num2str(z),' P=', num2str(P),'n=', num2str(n)); text (0, 0.1, actualtext);
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
        
        subplot (6,3,plotindex+1);
        cycbintimes=[repmat(sum(cycdur{un_indx}(:,1),1)./10, 1, th_nbins/2) repmat(sum(cycdur{un_indx}(:,2),1)./10, 1, th_nbins/2)];
        line ([binscale binscale+360], repmat(th_coupvect_LFPh(un_indx,:)./(cycbintimes),1,2));
        x=[0:1:720]; y=(cosd(x+180)+1)*(max(th_coupvect_LFPh(un_indx,:)./max(cycbintimes))).*0.5; line (x,y, 'Color', 'red');
        axis ([0 720 0 max(th_coupvect_LFPh(un_indx,:)./(th_time(un_indx)/20))*1.2]);
        set (gca, 'Xtick', 0:90:720); xlabel ('theta phase (deg)'); ylabel ('spike rate (s-1)')
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
        
        %phase plot for all spikes
        subplot (6,3,plotindex+2);
        plot (th_coupvect_LFP{un_indx}, [1:1:length(th_coupvect_LFP{un_indx})], 'b. ');
        axis ([0 360 0 length(th_coupvect_LFP{un_indx})]); xlabel ('phase (deg)'); ylabel ('spike #');
         %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
        
    end
end
clear th_nbins;

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
 
%saving the relevant data
datum = date;
writename=strcat(bs_name, bs_num, bs_exp, '_a_tht_cellcoupling');
save(writename, 'thcell','datum');
            
                     
        

        
           
            
           
           
 