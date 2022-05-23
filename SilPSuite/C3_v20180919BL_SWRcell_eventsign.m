%This code analyses SWR dependent firing: It does a shuffling for the SWR
%window to determine whether there are significantly more or less spikes
%found than expected in the surrounding time period

clear all; %close all

%USER INTERVENTION
bs_name='B'; bs_num='216'; bs_exp='b'; bs_typ='sp';           %unit basic data entry
fl_ord=[1];                                                 %file extension numbers belonging to the present cell
un_chans=[130:171];                       %channels containing the spike times for the units
un_num=[0, 4, 10, 11, 12, 20, 23, 47, 218, 35, 39, 40, 48, 49, 236, 54, 62, 190, 64, 65, 66, 82, 92, 103, 227, 104, 106, 115, 125, 126, 129, 132, 183, 215, 141, 144, 147, 149, 150, 152, 153, 154];                                               %unit numbers (in kwik file)
un_shnk =[repmat([1],1,9),repmat([2],1,6),repmat([3],1,3),repmat([4],1,3),repmat([5],1,4),repmat([6],1,3),repmat([7],1,6),repmat([8],1,8)];                  %tetrodes from which unit was isolated

%USER INTERVENTION
%set the numer of random placements to be performed for one SWR and the size of the random window to be included in the random placements (in sec).
Rnum = 10;
Rwin = 2;
%indicate one lfp channel for the file length extraction
ch_lfp=12;

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
end; clear fl_indx ch_lfp

%get the SWR periods together for the whole recording
sw_pers=[];
for fl_indx=1:length(fl_ord)
    fl_name=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'WAW.mat');   %file name
    sw_pname=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'SWP.txt');  %swr periods
    if exist(sw_pname,'file')==2;
        sw_persF=dlmread(sw_pname)+Fstarts(fl_indx); 
        if isempty(sw_pers)
            sw_pers=sw_persF;
            clear sw_persF
        else
            sw_pers=[sw_pers; sw_persF];
            clear sw_persF
        end
    end;
    clear sw_pname fl_name;
end; clear fl_indx;
sw_lengths=sw_pers(:,2)-sw_pers(:,1);
sw_cumlength=cumsum(sw_lengths);
sw_gaps=sw_pers(:,1)-[0; sw_cumlength(1:(end-1))];
%get the unit spike times together for the whole recording
spks{length(un_chans)}=[];
for un_indx=1:length(un_chans)
    for fl_indx=1:length(fl_ord); %file stepper
        fl_name=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'WAW.mat');   %file name
        %reading the spike sequence for the cell file
        ch_name=strcat(bs_name, bs_typ, bs_num, bs_exp,'_Ch', num2str(un_chans(un_indx)));
        ch_load=load(num2str(fl_name),num2str(ch_name));
        fl_spks=(ch_load.(ch_name).('times')); clear ch_name ch_load fl_name;
        fl_spks=fl_spks+Fstarts(fl_indx);
        if isempty(spks{un_indx})
            spks{un_indx}=fl_spks;
            clear fl_spks
        else
            spks{un_indx}=[spks{un_indx}; fl_spks];
            clear fl_spks
        end
    end; clear fl_indx
end;clear un_indx

%making of modified variables
for un_indx=1:length(un_chans)
    spks_all=spks{un_indx};
    spks_out=periodcut(spks_all, sw_pers(:,1:2), 0);
    spks_adj=zeros(size(spks_out));
    for je=1:size(sw_pers,1)
        spks_adj(spks_out>sw_pers(je,1))=sw_cumlength(je);
    end; clear je
    spks_rem=spks_out-spks_adj; clear spks_adj
    
    %extract data spk nums and shuffled numbers
    num_in{un_indx} = zeros(size(sw_pers,1),1);
    num_shf{un_indx}= zeros(size(sw_pers,1),Rnum);
    
    %stepper for swrs
    for je=1:size(sw_pers,1)
        %number of spikes in swr
        num_in{un_indx}(je)=length(spks_all(spks_all>sw_pers(je,1)&spks_all<=sw_pers(je,2)));
        
        %shuffling
        if sw_gaps(je,1)<Rwin/2
            for ka=1:Rnum
                sh_sta = (rand(1)*(sw_gaps(je)+(Rwin/2)-sw_lengths(je)));
                sh_end = sh_sta+sw_lengths(je);
                num_shf{un_indx}(je,ka)=length(spks_rem(spks_rem>sh_sta&spks_rem<=sh_end));
                clear sh_sta sh_end
            end; clear ka
        elseif sw_gaps(je,1)+sw_lengths(je,1)>(Fstarts(end)-sum(sw_lengths)-(Rwin/2))
            for ka=1:Rnum
                sh_sta = (Fstarts(end)-sum(sw_lengths))-(rand(1)*((Rwin/2)+(Fstarts(end)-sum(sw_lengths))-sw_gaps(je)));
                sh_end = sh_sta+sw_lengths(je);
                num_shf{un_indx}(je,ka)=length(spks_rem(spks_rem>sh_sta&spks_rem<=sh_end));
                clear sh_sta sh_end
            end; clear ka
        else
            for ka=1:Rnum
                sh_sta=(sw_gaps(je)-(rand(1)*(Rwin-sw_lengths(je))))+(Rwin/2);
                sh_end = sh_sta+sw_lengths(je);
                num_shf{un_indx}(je,ka)=length(spks_rem(spks_rem>sh_sta&spks_rem<=sh_end));
                clear sh_sta sh_end
            end; clear ka
        end
    end; clear je
end; clear un_indx


%----PLOTS
figindex=0;
for un_indx=1:length(un_chans)
    figindex_T=ceil(un_indx/6);
    if figindex_T~=figindex
        figindex=figindex+1;
        figname=strcat(bs_name, bs_num, bs_exp, num2str(figindex),'_SWRstats');
        figure ('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
        plotindex=1;
    else plotindex=((un_indx-((figindex-1)*6))*3)-2;
    end
    
    if (max(num_in{un_indx})+max(max(num_shf{un_indx})))>0
        %1st histogram - spike numbers in SWRs
        subplot (6,3,plotindex)
        Xaxis = [0:1:max([max(num_in{un_indx}) max(max(num_shf{un_indx}))])];
        hist(num_in{un_indx}, Xaxis);
        title ('distribution in SWRs'); xlabel ('spike count'); ylabel ('SWR count');
        axis ([-1 max([max(num_in{un_indx}) max(max(num_shf{un_indx}))])+1 0 1.1*(max(hist(num_in{un_indx}, Xaxis)))]);
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
        
        %2nd histogram - spike numbers in shuffling windows
        subplot(6,3,plotindex+1)
        Xaxis = [0:1:max([max(num_in{un_indx}) max(max(num_shf{un_indx}))])];
        bar(Xaxis,sum(hist(num_shf{un_indx}, Xaxis)'));
        title ('distribution in shuffled windows'); xlabel ('spike count'); ylabel ('window count');
        axis ([-1 max([max(num_in{un_indx}) max(max(num_shf{un_indx}))])+1 0 1.1*max(sum(hist(num_shf{un_indx}, Xaxis)'))]);
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
    end
    %textual results - statistic
    subplot(6,3,plotindex+2)
    axis off;
    actualtext= strcat('unit:', bs_name, bs_num, bs_exp, 't', num2str(un_shnk(un_indx)),'u',num2str(un_num(un_indx))); text(-0.2, 1, actualtext);
    actualtext= strcat('mean SWR duration:', num2str(mean(sw_lengths)*1000), 'ms'); text(-0.2, 0.88, actualtext);
    actualtext= strcat('in/out N:', num2str(sum(num_in{un_indx})),'/',num2str(mean(sum(num_shf{un_indx})))); text(-0.2, 0.76, actualtext);
    stati.data(1,un_indx)=sum(num_in{un_indx});
    stati.data(2,un_indx)=mean(sum(num_shf{un_indx}));
    actualtext= strcat('mean in SWR frq:', num2str(sum(num_in{un_indx})/sum(sw_lengths)), 'Hz'); text(-0.2, 0.64, actualtext);
    stati.data(3,un_indx)=(sum(num_in{un_indx})/sum(sw_lengths));
    actualtext= strcat('mean out SWR frq:', num2str((sum(sum(num_shf{un_indx})))/(sum(sw_lengths)*Rnum)), 'Hz'); text(-0.2, 0.52, actualtext);
    stati.data(4,un_indx)=((sum(sum(num_shf{un_indx})))/(sum(sw_lengths)*Rnum));
    [p,h,stats]=ranksum(num_in{un_indx}, reshape(num_shf{un_indx},size(num_shf{un_indx},1)*size(num_shf{un_indx},2),1));
    stati.data(6,un_indx)=p;
    actualtext=strcat('P/ranksum(Mann-Whitney U):', num2str(p),'/',num2str(stats.ranksum)); text(-0.2, 0.40, actualtext);
    [h,p]=kstest2(num_in{un_indx}, reshape(num_shf{un_indx},size(num_shf{un_indx},1)*size(num_shf{un_indx},2),1));
    stati.data(5,un_indx)=p;
    actualtext=strcat('P(Kolmogorov-Smirnov):', num2str(p)); text(-0.2, 0.28, actualtext);
    actualtext=strcat('shuffling number:', num2str(Rnum)); text(-0.2, 0.16, actualtext);
    actualtext=strcat('shuffling window:', num2str(Rwin), 's'); text(-0.2, 0.04, actualtext);
    %expand the graph
    pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
    set(gca, 'Clipping', 'off')
    stati.legend={'Nin','Nout','FRQin','FRQout','P(K-S)','P(M-W)'};
end

clear p h spks_all spks_out spks_rem sw_cumlength Xaxis actualtext figindex figindex_T figname plotindex 

datum=date;

%saving the relevant data
writename=strcat(bs_name, bs_num, bs_exp, '_a_swr_eventsignificance');
save(writename, 'Fstarts', 'Rnum','Rwin','num_in','num_shf','spks','sw_pers','sw_gaps','datum','stati');
 
clear writename figname
         
        
%2018.04.13 bug corrected in line 82:
%spks_rem=spks_out+spks_adj; clear spks_adj changed to spks_rem=spks_out-spks_adj; clear spks_adj
%2018.09.19 stati variable added to store data for easier recovery in summary excel files

       
       
        



            