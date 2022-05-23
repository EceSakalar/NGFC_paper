%This code analyses SWR dependent firing for cells recorded with the glass electrode 
%It does histograms and statistics for SWR dependent firing 
%and allows to include segments of the total recording

clear all; %close all

%USER INTERVENTION
%unit basic data entry
bs_name='B'; bs_num='183'; bs_exp='a'; bs_typ='sx';

%unit firing information content
fl_ord=[1];                                       %file extension numbers belonging to the present cell
gl_chan=[22];                                     %channel containing the spike times for the cell
gl_num=[1];                                       %cell number to be analysed numbers (in kwik file)
gl_limS=[2682];                                   %recording of cell starts (in seconds, for each file)
gl_limE=[2931];                                   %recording of cell starts (in seconds, for each file)


%analysis information HISTOGRAM
%give the number of bins for the center (middle) and surroundings (side) of SWRs histogram
middle=8; side=12;
%give the binsize (ms) and the time before and after (sammple number) for the startaligned and endaligned histograms
binsize=20; binsize=binsize/1000; before=20; after=20;

%analysis information STATISTICS
%set the numer of random placements to be performed for one SWR and the size of the random window to be included in the random placements (in sec).
Rnum = 10;
Rwin = 2;
%indicate one lfp channel for the file length extraction
ch_lfp=15;

%----------HISTOGRAM SCRIPT STARTS
SWR_spikecount=[];

for fl_indx=1:length(fl_ord);                                              %file stepper
    fl_name=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'WAW.mat');   %file name
    sw_pname=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'SWP.txt');  %swr periods
    
    %reading the spike sequence for the cell
    ch_name=strcat(bs_name, bs_typ, bs_num, bs_exp,'_Ch', num2str(gl_chan));
    ch_load=load(num2str(fl_name),num2str(ch_name));
    fl_spks=(ch_load.(ch_name).('times')); clear ch_name ch_load fl_name;
        
    if exist(sw_pname,'file')==2; sw_pers=dlmread(sw_pname); clear sw_pname;
        sw_pers(sw_pers(:,1)<gl_limS(fl_indx),:)=[];
        sw_pers(sw_pers(:,2)>gl_limE(fl_indx),:)=[];
        
        fl_swrspks=periodcut(fl_spks,sw_pers(:,1:2),1);
        
        %I calculate the SWR time spent in this file
        SWR_totaltime_file=sum(((sw_pers(:,2))-(sw_pers(:,1))),1);
        if fl_indx==1
            SWR_totaltime=SWR_totaltime_file;
        else
            SWR_totaltime=SWR_totaltime+SWR_totaltime_file; 
        end; clear SWR_totaltime_file;
            
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
        if isempty(SWR_spikecount)
            SWR_startaligned=startaligned_f;
            SWR_endaligned=endaligned_f;
            SWR_spikecount=SWR_spikecount_f;
            SWR_bintime=SWR_bintime_f;
            SWR_spikecounttoblank=SWR_spikecounttoblank_f;
            SWR_bintimetoblank=SWR_bintimetoblank_f;
        else
            SWR_startaligned=[SWR_startaligned; startaligned_f];
            SWR_endaligned=[SWR_endaligned; endaligned_f];
            SWR_spikecount=[SWR_spikecount; SWR_spikecount_f];
            SWR_bintime=[SWR_bintime; SWR_bintime_f];
            SWR_spikecounttoblank=[SWR_spikecounttoblank; SWR_spikecounttoblank_f];
            SWR_bintimetoblank=[SWR_bintimetoblank; SWR_bintimetoblank_f];
        end
        clear SWR_spikecount_f SWR_bintime_f SWR_spikecounttoblank_f SWR_bintimetoblank_f SWR_binstarts_f SWR_binends_f
        clear sw_pers fl_swrspks startaligned_f endaligned_f
    end; clear fl_spks
end; clear fl_indx

SWR_spikecountblanked=SWR_spikecount-SWR_spikecounttoblank;
SWR_bintimeblanked=SWR_bintime-SWR_bintimetoblank;
clear SWR_bintimetoblank SWR_spikecounttoblank


SWR.Xaxis=[1:1:middle+(2*side)];
SWR.totaltime=SWR_totaltime;
SWR.aligned_xaxis=[((-1*(binsize*(before-0.5)))*1000):(binsize*1000):(((after-0.5)*binsize)*1000)];

%plotting the results

figname=strcat(bs_name, bs_num, bs_exp, 'gl', num2str(gl_num),'_autocorrelogram');
figure ('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')

subplot (2,3,1);
%textual results
if max(max(SWR_spikecountblanked))>0
    axis off
    text(0,1,strcat(bs_name, bs_num, bs_exp, ' glass', num2str(gl_num)))
    %calculus
    SWR.spikenumberH=sum(SWR_spikecountblanked,1);
    SWR.number=size(SWR_spikecount,1);
    SWR.activenumber=size(find((sum((SWR_spikecountblanked(:,(side+1):(side+middle))),2))>0),1);
    SWR.spikenumber=sum(SWR.spikenumberH((side+1):(side+middle)));
    SWR.spikefrequency=SWR.spikenumber/SWR.totaltime;
    %number of SWRs
    actualtext= strcat('Number of SWRs analysed: ', num2str (SWR.number)); text (0, 0.8, actualtext);
    %number of active SWRs
    actualtext= strcat('Number of active SWRs : ', num2str (SWR.activenumber)); text (0, 0.6, actualtext);
    %mean number of APs in SWRs
    actualtext= strcat('Mean number of APs in an SWR:', num2str((SWR.spikenumber)/(SWR.activenumber))); text (0, 0.4, actualtext);
    %mean AP frequency
    actualtext=strcat('Mean AP frequency under SWR:', num2str (SWR.spikefrequency)); text (0, 0.2, actualtext);
    %overall time spent in SWR
    actualtext=strcat('overall SWR time:', num2str (SWR.totaltime), ' s'); text (0, 0, actualtext);
    
    %SWR spike frq histrogram
    subplot(2,3,2)
    SWR.spikefrqH=(sum(SWR_spikecountblanked,1))./(sum(SWR_bintimeblanked,1));
    bar (SWR.Xaxis, SWR.spikefrqH);
    set (gca, 'Xtick', 0:1:middle+(2*side));
    axis ([0 (middle+(2*side)+1) 0 max(SWR.spikefrqH)*1.1]);
    xlabel('phase (of SWR)'); ylabel('spike/s');
    line([side+0.5 side+0.5],[0 max(SWR.spikefrqH)*1.1], 'Linewidth',2,'Color','red');
    line([side+middle+0.5 side+middle+0.5],[0 max(SWR.spikefrqH)*1.1], 'Linewidth',2,'Color','red');
       
    %histograms aligned to the start or to the end of SWRs
    subplot(2,3,3)
    SWR_startaligned(:,1)=[];SWR_startaligned(:,end)=[];
    SWR_endaligned(:,1)=[]; SWR_endaligned(:,end)=[];
    SWR.startaligned=sum(SWR_startaligned,1)./(binsize*size(SWR_startaligned,1));
    SWR.endaligned=sum(SWR_endaligned,1)./(binsize*size(SWR_endaligned,1));
    line(SWR.aligned_xaxis,SWR.startaligned,'Linewidth',2,'Color','red');line(SWR.aligned_xaxis,SWR.endaligned,'Linewidth',2,'Color','blue');
    set(gca, 'Xtick', (before*binsize*-1000):binsize*2000:(after*binsize*1000));
    axis([min(SWR.aligned_xaxis) max(SWR.aligned_xaxis) 0 (max([max(SWR.startaligned) max(SWR.endaligned)]))*1.1]);
    xlabel ('time (ms)'); ylabel ('spike rate(s-1)');  grid on;
end
clear SWR_startaligned SWR_endaligned
clear SWR_spikecountblanked SWR_bintimeblanked SWR_spikecounttoblank SWR_bintimetoblank SWR_bintime SWR_spikecount
clear after before binsize middle side actualtext
clear SWR_totaltime

%----------STATISTICS SCRIPT STARTS
%get the file lengths (starting points in seconds in one vector)
Fstarts=0;
for fl_indx=1:length(fl_ord)
    fl_name=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'WAW.mat');   %file name
    ch_name=strcat(bs_name, bs_typ, bs_num, bs_exp,'_Ch', num2str(ch_lfp));
    ch_load=load(num2str(fl_name),num2str(ch_name));
    Fdata=(ch_load.(ch_name).('length'));
    Finter=(ch_load.(ch_name).('interval'));
    Flength=Fdata*Finter; clear Fdata Finter h_load ch_name fl_name;
    Fstarts=[Fstarts Fstarts(fl_indx)+Flength]; clear Flength;
end; clear fl_indx ch_lfp ch_load

%get the SWR periods together for the whole recording
sw_pers=[];
for fl_indx=1:length(fl_ord)
    sw_pname=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'SWP.txt');  %swr periods
    if exist(sw_pname,'file')==2;
        sw_persF=dlmread(sw_pname)+Fstarts(fl_indx);
        sw_persF(sw_persF(:,1)<(gl_limS(fl_indx)+Fstarts(fl_indx)),:)=[];
        sw_persF(sw_persF(:,2)>(gl_limE(fl_indx)+Fstarts(fl_indx)),:)=[];
        
        if isempty(sw_pers)
            sw_pers=sw_persF;
        else
            sw_pers=[sw_pers; sw_persF];
        end; clear sw_persF
    end;
    clear sw_pname fl_name;
end; clear fl_indx;
sw_lengths=sw_pers(:,2)-sw_pers(:,1);
sw_cumlength=cumsum(sw_lengths);
sw_gaps=sw_pers(:,1)-[0; sw_cumlength(1:(end-1))];
%get the unit spike times together for the whole recording
spks=[];
for fl_indx=1:length(fl_ord); %file stepper
    fl_name=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'WAW.mat');   %file name
    %reading the spike sequence for the cell file
    ch_name=strcat(bs_name, bs_typ, bs_num, bs_exp,'_Ch', num2str(gl_chan));
    ch_load=load(num2str(fl_name),num2str(ch_name));
    fl_spks=(ch_load.(ch_name).('times')); clear ch_name ch_load fl_name;
    fl_spks=fl_spks+Fstarts(fl_indx);
    if isempty(spks)
        spks=fl_spks;
    else
        spks=[spks; fl_spks];
    end; clear fl_spks
end; clear fl_indx

%making of modified variables
spks_all=spks;
spks_out=periodcut(spks_all, sw_pers(:,1:2), 0);
spks_adj=zeros(size(spks_out));
for je=1:size(sw_pers,1)
    spks_adj(spks_out>sw_pers(je,1))=sw_cumlength(je);
end; clear je
spks_rem=spks_out-spks_adj; clear spks_adj

%extract data spk nums and shuffled numbers
num_in = zeros(size(sw_pers,1),1);
num_shf= zeros(size(sw_pers,1),Rnum);

%stepper for swrs
for je=1:size(sw_pers,1)
    %number of spikes in swr
    num_in(je)=length(spks_all(spks_all>sw_pers(je,1)&spks_all<=sw_pers(je,2)));
    
    %shuffling
    if sw_gaps(je,1)<Rwin/2
        for ka=1:Rnum
            sh_sta = (rand(1)*(sw_gaps(je)+(Rwin/2)-sw_lengths(je)));
            sh_end = sh_sta+sw_lengths(je);
            num_shf(je,ka)=length(spks_rem(spks_rem>sh_sta&spks_rem<=sh_end));
            clear sh_sta sh_end
        end; clear ka
    elseif sw_gaps(je,1)+sw_lengths(je,1)>(Fstarts(end)-sum(sw_lengths)-(Rwin/2))
        for ka=1:Rnum
            sh_sta = (Fstarts(end)-sum(sw_lengths))-(rand(1)*((Rwin/2)+(Fstarts(end)-sum(sw_lengths))-sw_gaps(je)));
            sh_end = sh_sta+sw_lengths(je);
            num_shf(je,ka)=length(spks_rem(spks_rem>sh_sta&spks_rem<=sh_end));
            clear sh_sta sh_end
        end; clear ka
    else
        for ka=1:Rnum
            sh_sta=(sw_gaps(je)-(rand(1)*(Rwin-sw_lengths(je))))+(Rwin/2);
            sh_end = sh_sta+sw_lengths(je);
            num_shf(je,ka)=length(spks_rem(spks_rem>sh_sta&spks_rem<=sh_end));
            clear sh_sta sh_end
        end; clear ka
    end
end; clear je

%----PLOTS
if (max(num_in)+max(max(num_shf)))>0
    %1st histogram - spike numbers in SWRs
    subplot (2,3,5)
    Xaxis = [0:1:max([max(num_in) max(max(num_shf))])];
    hist(num_in, Xaxis);
    title ('distribution in SWRs'); xlabel ('spike count'); ylabel ('SWR count');
    axis ([-1 max([max(num_in) max(max(num_shf))])+1 0 1.1*(max(hist(num_in, Xaxis)))]);
    
    %2nd histogram - spike numbers in shuffling windows
    subplot(2,3,6)
    Xaxis = [0:1:max([max(num_in) max(max(num_shf))])];
    bar(Xaxis,sum(hist(num_shf, Xaxis)'));
    title ('distribution in shuffled windows'); xlabel ('spike count'); ylabel ('window count');
    axis ([-1 max([max(num_in) max(max(num_shf))])+1 0 1.1*max(sum(hist(num_shf, Xaxis)'))]);
end
%textual results - statistic
subplot(2,3,4)
axis off;
actualtext= strcat('unit:', bs_name, bs_num, bs_exp, 'glass', num2str(gl_num)); text(-0.2, 1, actualtext);
actualtext= strcat('mean SWR duration:', num2str(mean(sw_lengths)*1000), 'ms'); text(-0.2, 0.88, actualtext);
actualtext= strcat('in/out N:', num2str(sum(num_in)),'/',num2str(mean(sum(num_shf)))); text(-0.2, 0.76, actualtext);
actualtext= strcat('mean in SWR frq:', num2str(sum(num_in)/sum(sw_lengths)), 'Hz'); text(-0.2, 0.64, actualtext);
actualtext= strcat('mean out SWR frq:', num2str((sum(sum(num_shf)))/(sum(sw_lengths)*Rnum)), 'Hz'); text(-0.2, 0.52, actualtext);
[p,h,stati]=ranksum(num_in, reshape(num_shf,size(num_shf,1)*size(num_shf,2),1));
actualtext=strcat('P/ranksum(Mann-Whitney U):', num2str(p),'/',num2str(stati.ranksum)); text(-0.2, 0.40, actualtext);
[h,p]=kstest2(num_in, reshape(num_shf,size(num_shf,1)*size(num_shf,2),1));
actualtext=strcat('P(Kolmogorov-Smirnov):', num2str(p)); text(-0.2, 0.28, actualtext);
actualtext=strcat('shuffling number:', num2str(Rnum)); text(-0.2, 0.16, actualtext);
actualtext=strcat('shuffling window:', num2str(Rwin), 's'); text(-0.2, 0.04, actualtext);

clear p h spks_all spks_out spks_rem sw_cumlength Xaxis actualtext figindex figindex_T figname plotindex stati sw_length

%-----SAVE DATA
SWRstat.Fstarts=Fstarts;
SWRstat.Rnum=Rnum;
SWRstat.Rwin=Rwin;
SWRstat.numIN=num_in;
SWRstat.numSHF=num_shf;
SWRstat.spks=spks;
SWRstat.swPERS=sw_pers;
SWRstat.swGAPS=sw_gaps;

clear Fstarts Rnum Rwin num_in num_shf spks sw_pers sw_gaps

datum=date;

%saving the relevant data
writename=strcat(bs_name, bs_num, bs_exp, '_a_swr_GLASSUcoupling_u',num2str(gl_num));
save(writename, 'SWR', 'SWRstat','datum');

clear writename figname
        
        
       
       
        



            