%This code analyses SWR dependent firing: It does histograms for SWR dependent firing 

clear all; %close all

%USER INTERVENTION
bs_name='B'; bs_num='216'; bs_exp='b'; bs_typ='sp'; %unit basic data entry

%unit firing information content
fl_ord=[1];                                       %file extension numbers belonging to the present cell

un_chans =[130:171];                                                           %channels containing the spike times for the units
un_nums = [0, 4, 10, 11, 12, 20, 23, 47, 218, 35, 39, 40, 48, 49, 236, 54, 62, 190, 64, 65, 66, 82, 92, 103, 227, 104, 106, 115, 125, 126, 129, 132, 183, 215, 141, 144, 147, 149, 150, 152, 153, 154];          %unit numbers (in kwik file)
%shanks from which unit was isolated; shank 0 is glass electrode and allows shortening of the recording time
un_shnk = [repmat([1],1,9),repmat([2],1,6),repmat([3],1,3),repmat([4],1,3),repmat([5],1,4),repmat([6],1,3),repmat([7],1,6),repmat([8],1,8)];                       

%analysis information
%give the number of bins for the center (middle) and surroundings (side) of SWRs histogram
middle=8; side=12;
%give the binsize (ms) and the time before and after (sammple number) for the startaligned and endaligned histograms
binsize=20; binsize=binsize/1000; before=20; after=20;

%script starts
SWR_totaltime=0;
    
for un_indx=1:length(un_chans)
    SWR_spikecount{un_indx}=[];
end
for un_indx=1:length(un_chans)
    
    for fl_indx=1:length(fl_ord); %file stepper
        fl_name=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'WAW.mat');   %file name
        sw_pname=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'SWP.txt');  %swr periods
        
        %reading the spike sequence for the cell
        ch_name=strcat(bs_name, bs_typ, bs_num, bs_exp,'_Ch', num2str(un_chans(un_indx)));
        ch_load=load(num2str(fl_name),num2str(ch_name));
        fl_spks=(ch_load.(ch_name).('times')); clear ch_name ch_load fl_name;
        
        if exist(sw_pname,'file')==2; sw_pers=dlmread(sw_pname); clear sw_pname;
            
            fl_swrspks=periodcut(fl_spks,sw_pers(:,1:2),1);
            
            %I calculate the SWR time spent in this file
            if un_indx==1
                SWR_totaltime_file=sum(((sw_pers(:,2))-(sw_pers(:,1))),1);
                SWR_totaltime=SWR_totaltime+SWR_totaltime_file; clear SWR_totaltime_file;
            end
            
            %I calculate the bin sizes
            bin_start=(sw_pers(:,3)-sw_pers(:,1))/(middle/2); bin_end=(sw_pers(:,2)-sw_pers(:,3))/(middle/2); bin_side=(bin_start+bin_end)/2;
            
            %I calculate binmatrix for the four periods (before the SWR, first and second SWR halves and after the SWR)
            SWR_binstarts_k=repmat((sw_pers(:,1)-(bin_side*side)),1,(middle+2*side));
            SWR_binstarts_a=cumsum([zeros(size(bin_side,1),1) repmat(bin_side,1,side) repmat(bin_start,1,middle/2) repmat(bin_end,1,middle/2) repmat(bin_side,1,(side-1))],2);
            SWR_binends_a=cumsum([repmat(bin_side,1,side) repmat(bin_start,1,middle/2) repmat(bin_end,1,middle/2) repmat(bin_side,1,side)],2);
            SWR_binstarts_f=SWR_binstarts_k+SWR_binstarts_a;
            SWR_binends_f=SWR_binstarts_k+SWR_binends_a;
            clear bin_start bin_end bin_side SWR_binstarts_k SWR_binstarts_a SWR_binends_a;
            
            %I calculate the spike numbers for all the bins
            for j=1:size(SWR_binstarts_f,1)
                for k=1:size(SWR_binstarts_f,2)
                    SWR_spikecount_f(j,k)=length(find(fl_spks>=SWR_binstarts_f(j,k)&fl_spks<SWR_binends_f(j,k)));
                end
            end; clear j k
            SWR_bintime_f=SWR_binends_f-SWR_binstarts_f;
            
            SWR_spikecounttoblank_f=zeros(size(SWR_spikecount_f));
            SWR_bintimetoblank_f=zeros(size(SWR_bintime_f));
            
            %I'm blanking the spikes included in neighbouring SWRs
            for j=1:size(SWR_binstarts_f,1)
                for k=1:side
                    SWR_spikecounttoblank_f(j,k)=length(find(fl_swrspks>=SWR_binstarts_f(j,k)&fl_swrspks<SWR_binends_f(j,k)));
                end; clear k
                for k=(side+middle+1):(side+middle+side)
                    SWR_spikecounttoblank_f(j,k)=length(find(fl_swrspks>SWR_binstarts_f(j,k)&fl_swrspks<=SWR_binends_f(j,k)));
                end; clear k
            end; clear j
            
            %I'm blanking the time included in neighbouring SWRs
            for j=1:size(SWR_binstarts_f,1)
                if j>1
                    for k=1:side
                        ind=find(sw_pers(1:j-1,1)<SWR_binends_f(j,k)&sw_pers(1:j-1,2)>SWR_binstarts_f(j,k));
                        if size(ind,1)>0
                            for ind2=1:size(ind,1)
                                if SWR_binstarts_f(j,k)>sw_pers(ind(ind2),1)&&SWR_binends_f(j,k)<sw_pers(ind(ind2),2);
                                    SWR_bintimetoblank_f(j,k)=SWR_binends_f(j,k)-SWR_binstarts_f(j,k);
                                elseif SWR_binstarts_f(j,k)<=sw_pers(ind(ind2),1)
                                    SWR_bintimetoblank_f(j,k)=SWR_bintimetoblank_f(j,k)+(SWR_binends_f(j,k)-sw_pers(ind(ind2),1));
                                elseif SWR_binends_f(j,k)>=sw_pers(ind(ind2),2)
                                    SWR_bintimetoblank_f(j,k)=SWR_bintimetoblank_f(j,k)+(sw_pers(ind(ind2),2)-SWR_binstarts_f(j,k));
                                end
                            end; clear ind2;
                        end
                    end; clear ind k;
                end
                
                if j<size(SWR_binstarts_f, 1);
                    for k=(middle+side+1):(middle+(2*side))
                        ind=find(sw_pers(:,1)<SWR_binends_f(j,k)&sw_pers(:,2)>SWR_binstarts_f(j,k));
                        if size(ind,1)>0
                            for ind2=1:size(ind,1)
                                if SWR_binstarts_f(j,k)>sw_pers(ind(ind2),1)&&SWR_binends_f(j,k)<sw_pers(ind(ind2),2);
                                    SWR_bintimetoblank_f(j,k)=SWR_binends_f(j,k)-SWR_binstarts_f(j,k);
                                elseif SWR_binstarts_f(j,k)<=sw_pers(ind(ind2),1)
                                    SWR_bintimetoblank_f(j,k)=SWR_bintimetoblank_f(j,k)+(SWR_binends_f(j,k)-sw_pers(ind(ind2),1));
                                elseif SWR_binends_f(j,k)>=sw_pers(ind(ind2),2)
                                    SWR_bintimetoblank_f(j,k)=SWR_bintimetoblank_f(j,k)+(sw_pers(ind(ind2),2)-SWR_binstarts_f(j,k));
                                end
                            end; clear ind2
                        end
                    end; clear ind k
                end
            end; clear j;
            
            %calculating the SWR-instant spike-rate plots
            for j=1:size(sw_pers,1)
                sw_actspks=periodcut(fl_spks,sw_pers(j,:),1);
                if length(sw_actspks)>1
                    sw_actspks_pos=(sw_actspks-sw_pers(j,1))./(sw_pers(j,2)-sw_pers(j,1));
                    for k=1:length(sw_actspks)-1
                        sw_ratedat_per(k,1)=(sw_actspks_pos(k)+sw_actspks_pos(k+1))/2;
                        sw_ratedat_per(k,2)=1/(sw_actspks(k+1)-sw_actspks(k));
                    end; clear sw_actspks_pos k
                    
                    if exist('sw_ratedat','var')
                        sw_ratedat=[sw_ratedat; sw_ratedat_per];
                    else
                        sw_ratedat=sw_ratedat_per;
                    end; clear sw_ratedat_per;
                end; clear sw_actspks;
            end; clear j;
            
            %generating the endaligned and the startaligned binned histograms
            aligned_binner=repmat((-1*before*binsize)-(binsize/2):binsize:(after*binsize)+(binsize/2),size(sw_pers,1),1);
            start_binner=repmat(sw_pers(:,1),1,size(aligned_binner,2));
            end_binner=repmat(sw_pers(:,2),1,size(aligned_binner,2));
            startaligned_binner=start_binner+aligned_binner;
            endaligned_binner=end_binner+aligned_binner;
            startaligned_f=zeros(size(startaligned_binner)); endaligned_f=zeros(size(endaligned_binner));
            clear aligned_binner start_binner end_binner
            %actual binning
            for j=1:size(sw_pers,1)
                startaligned_f(j,:)=hist(fl_spks,startaligned_binner(j,:));
                endaligned_f(j,:)=hist(fl_spks,endaligned_binner(j,:));
            end; clear j startaligned_binner endaligned_binner
            
            %I'm gathering the data for the actual unit SWR
            if isempty(SWR_spikecount{un_indx})
                SWR_startaligned{un_indx}=startaligned_f;
                SWR_endaligned{un_indx}=endaligned_f;
                SWR_spikecount{un_indx}=SWR_spikecount_f;
                SWR_bintime{un_indx}=SWR_bintime_f;
                SWR_spikecounttoblank=SWR_spikecounttoblank_f;
                SWR_bintimetoblank=SWR_bintimetoblank_f;
            else
                SWR_startaligned{un_indx}=[SWR_startaligned{un_indx}; startaligned_f];
                SWR_endaligned{un_indx}=[SWR_endaligned{un_indx}; endaligned_f];
                SWR_spikecount{un_indx}=[SWR_spikecount{un_indx}; SWR_spikecount_f];
                SWR_bintime{un_indx}=[SWR_bintime{un_indx}; SWR_bintime_f];
                SWR_spikecounttoblank=[SWR_spikecounttoblank; SWR_spikecounttoblank_f];
                SWR_bintimetoblank=[SWR_bintimetoblank; SWR_bintimetoblank_f];
            end
            clear SWR_spikecount_f SWR_bintime_f SWR_spikecounttoblank_f SWR_bintimetoblank_f SWR_binstarts_f SWR_binends_f
            clear sw_pers fl_swrspks startaligned_f endaligned_f
        end; clear fl_spks
    end; clear fl_indx
    
    SWR_spikecountblanked{un_indx}=SWR_spikecount{un_indx}-SWR_spikecounttoblank;
    SWR_bintimeblanked{un_indx}=SWR_bintime{un_indx}-SWR_bintimetoblank;
    clear SWR_bintimetoblank SWR_spikecounttoblank
end

SWR.Xaxis=[1:1:middle+(2*side)];
SWR.totaltime=SWR_totaltime;
SWR.aligned_xaxis=[((-1*(binsize*(before-0.5)))*1000):(binsize*1000):(((after-0.5)*binsize)*1000)];
%plotting the results
figindex=0;
for un_indx=1:length(un_chans)
    figindex_T=ceil(un_indx/6);
    if figindex_T~=figindex
        figindex=figindex+1;
        figname=strcat(bs_name, bs_num, bs_exp, num2str(figindex),'_autocorrelogram');
        figure ('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
        plotindex=1;
    else plotindex=((un_indx-((figindex-1)*6))*3)-2;
    end
    subplot (6,3,plotindex);
    %textual results
    if max(max(SWR_spikecountblanked{un_indx}))>0
        axis off
        text(0,1,strcat(bs_name, bs_num, bs_exp, ' t', num2str(un_shnk(un_indx)),'u',num2str(un_nums(un_indx))))
        %calculus
        SWR.spikenumberH{un_indx}=sum(SWR_spikecountblanked{un_indx},1);
        SWR.number{un_indx}=size(SWR_spikecount{un_indx},1);
        SWR.activenumber{un_indx}=size(find((sum((SWR_spikecountblanked{un_indx}(:,(side+1):(side+middle))),2))>0),1);
        SWR.spikenumber{un_indx}=sum(SWR.spikenumberH{un_indx}((side+1):(side+middle)));
        SWR.spikefrequency{un_indx}=SWR.spikenumber{un_indx}/SWR.totaltime;
        %number of SWRs
        actualtext= strcat('Number of SWRs analysed: ', num2str (SWR.number{un_indx})); text (0, 0.8, actualtext);
        %number of active SWRs
        actualtext= strcat('Number of active SWRs : ', num2str (SWR.activenumber{un_indx})); text (0, 0.6, actualtext);
        %mean number of APs in SWRs
        actualtext= strcat('Mean number of APs in an SWR:', num2str((SWR.spikenumber{un_indx})/(SWR.activenumber{un_indx}))); text (0, 0.4, actualtext);
        %mean AP frequency
        actualtext=strcat('Mean AP frequency under SWR:', num2str (SWR.spikefrequency{un_indx})); text (0, 0.2, actualtext);
        %overall time spent in SWR
        actualtext=strcat('overall SWR time:', num2str (SWR.totaltime), ' s'); text (0, 0, actualtext);
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
        
        %SWR spike frq histrogram
        subplot(6,3,plotindex+1)
        SWR.spikefrqH{un_indx}=(sum(SWR_spikecountblanked{un_indx},1))./(sum(SWR_bintimeblanked{un_indx},1));
        bar (SWR.Xaxis, SWR.spikefrqH{un_indx});
        set (gca, 'Xtick', 0:1:middle+(2*side));
        axis ([0 (middle+(2*side)+1) 0 max(SWR.spikefrqH{un_indx})*1.1]);
        xlabel('phase (of SWR)'); ylabel('spike/s');
        line([side+0.5 side+0.5],[0 max(SWR.spikefrqH{un_indx})*1.1], 'Linewidth',2,'Color','red');
        line([side+middle+0.5 side+middle+0.5],[0 max(SWR.spikefrqH{un_indx})*1.1], 'Linewidth',2,'Color','red');
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
        
        %histograms aligned to the start or to the end of SWRs
        subplot(6,3,plotindex+2)
        SWR_startaligned{un_indx}(:,1)=[];SWR_startaligned{un_indx}(:,end)=[];
        SWR_endaligned{un_indx}(:,1)=[]; SWR_endaligned{un_indx}(:,end)=[];
        SWR.startaligned{un_indx}=sum(SWR_startaligned{un_indx},1)./(binsize*size(SWR_startaligned{un_indx},1));
        SWR.endaligned{un_indx}=sum(SWR_endaligned{un_indx},1)./(binsize*size(SWR_endaligned{un_indx},1));
        line(SWR.aligned_xaxis,SWR.startaligned{un_indx},'Linewidth',2,'Color','red');line(SWR.aligned_xaxis,SWR.endaligned{un_indx},'Linewidth',2,'Color','blue');
        set(gca, 'Xtick', (before*binsize*-1000):binsize*2000:(after*binsize*1000));
        axis([min(SWR.aligned_xaxis) max(SWR.aligned_xaxis) 0 (max([max(SWR.startaligned{un_indx}) max(SWR.endaligned{un_indx})]))*1.1]);
        xlabel ('time (ms)'); ylabel ('spike rate(s-1)');  grid on;
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
    end
end

clear SWR_startaligned SWR_endaligned
clear SWR_spikecountblanked SWR_bintimeblanked SWR_spikecounttoblank SWR_bintimetoblank SWR_bintime SWR_spikecount
clear after before binsize middle side actualtext

datum=date;

%saving the relevant data
writename=strcat(bs_name, bs_num, bs_exp, '_a_swr_eventcoupling');
save(writename, 'SWR', 'datum');

clear writename figname
        
        
       
       
        



            