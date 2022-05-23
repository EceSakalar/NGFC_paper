%This program analyses the gamma oscillatory coupling relative to CSD
%traces (spectral aproach; multichannel LFP and CSD traces) 
clear all; close all

bs_name='B'; bs_num='178'; bs_exp='g'; bs_typ='sx';                 %unit basic data entry

ch_n=16;                                                            %number of channels
ch_rate=2000;                                                       %give the silicon probe sampling rate (in Hz)

un_chans=[19:38 42:44];                                                   %channels containing the spike times for the units
un_nums=[5 11 26 33 52 64 69 73 75 83 85 88 89 96 97 98 100 103 104 106 109 110 1];                 %unit numbers (in kwik file)
un_shank =[repmat([],1,12) repmat([2],1,22) 0];                        %shanks from which unit was isolated

ch_def=[];                                                          %defective channel (one is standing for the ventralmost, and the linear number of channel is required not the identifier
fl_ord=[1];                                                         %file extension numbers belonging to the present cell

%gamma wavelet properties
gm_frsta=250; gm_frfin=20; gm_cfnum=69; gm_wavelet{1}='cmor1-1.5'; gm_plimit=0.05;
gm_frcnt=centfrq(gm_wavelet{1}); gm_scsta=gm_frcnt/(gm_frsta/ch_rate); gm_scfin=gm_frcnt/(gm_frfin/ch_rate); gm_scint=(gm_scfin-gm_scsta)/gm_cfnum; gm_fraxis=gm_frcnt./((gm_scsta:gm_scint:gm_scfin)./ch_rate);
gm_nbins=20;gm_binsize=360/gm_nbins;

for je=1:ch_n
    for un_indx=1:length(un_chans)
        gm_cellcwtspect{je,un_indx}=[];
    end; clear un_indx
end; clear je

for je=1:ch_n;                                                      %CSD channel stepper
    for fl_indx=1:length(fl_ord);                                   %file stepper
        fl_name=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'WAW.mat'); %file name
        th_pname=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'THP.txt');%theta periods
        
        if exist(th_pname,'file')==2; th_pers=dlmread(th_pname); clear th_pname;
            
            %load the csd traces
            csd_name=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'CSD.mat'); load(csd_name);
            
            %check if it is edge and if it is defected channel
            if ~isnan(CSD(je,1))&&isempty(ch_def(ch_def==(je)))
                %invert the CSD so that the sink is down and the source is up
                CSDchan=CSD(je,:)*-1; clear CSD csd_name;
                
                %reading the spike sequence for the cells, transforming it to sample number instead of time
                for un_indx=1:length(un_chans)
                    ch_name=strcat(bs_name, bs_typ, bs_num, bs_exp,'_Ch', num2str(un_chans(un_indx)));
                    ch_load=load(num2str(fl_name),num2str(ch_name));
                    fl_spks{un_indx}=((ch_load.(ch_name).('times')).*ch_rate); clear ch_name ch_load;
                end
                
                %calculate the theta segment wavelet transforms
                for ka=1:size(th_pers,1); %theta period stepper
                    sg_start=(th_pers(ka,1)*ch_rate); sg_end=(th_pers(ka,2)*ch_rate);
                    sg_CSD=CSDchan((round(sg_start)-0.5*ch_rate):(round(sg_end)+0.5*ch_rate));
                    %calculating the wavelet transform for the actual channel of the actual segment.
                    sg_cwt=cwt(sg_CSD, gm_scsta:gm_scint:gm_scfin, gm_wavelet{1});
                    sg_cwt=conj(sg_cwt);
                    sg_cwt=sg_cwt(:,(1+0.5*ch_rate):(1+0.5*ch_rate+round(sg_end)-round(sg_start)));
                    
                    %extract spike triggered CWT for the particular segment for all units and add to a 2D cell array
                    for un_indx=1:length(un_chans)
                        sg_spks=round((fl_spks{un_indx}(fl_spks{un_indx}>sg_start&fl_spks{un_indx}<sg_end)-sg_start+1)');
                        sg_gspect=sg_cwt(:,sg_spks);
                        if isempty(gm_cellcwtspect{je,un_indx})
                            gm_cellcwtspect{je,un_indx}=sg_gspect;
                        else
                            gm_cellcwtspect{je,un_indx}=[gm_cellcwtspect{je,un_indx} sg_gspect];
                        end;
                        clear sg_spks sg_gspect
                    end; clear un_indx
                    clear sg_cwt sg_start sg_end sg_CSD
                end; clear ka th_pers CSDchan
                strcat('CSD_channel_', num2str(je), '_done for all units')
                clear fl_spks
            else
                strcat('CSD_channel_', num2str(je), '_is_not included')
            end
        end
    end; clear fl_indx
end; clear je
    




%%%%%%%%% ITT TARTOK
    


%plotting results gamma coupling (during theta)
%calculus
    for un_indx=1:length(un_chans)
        gm_cellcwtspect_r{un_indx}(1:ch_n,1:length(gm_fraxis))=NaN;
        gm_cellcwtspect_m{un_indx}(1:ch_n,1:length(gm_fraxis))=NaN;
        gm_cellcwtspect_p{un_indx}(1:ch_n,1:length(gm_fraxis))=NaN;
        gm_cellcwtspect_z{un_indx}(1:ch_n,1:length(gm_fraxis))=NaN;
        gm_cellcwtspect_sign{un_indx}(1:ch_n,1:length(gm_fraxis))=NaN;
    end

for un_indx=1:length(un_chans)
    %histograms
    figname=strcat(bs_name, bs_num, bs_exp, ' shk', num2str(un_shank(un_indx)),'u',num2str(un_nums(un_indx)));
    figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
    subplot(15,3,1);
    axis off;
    text(0,1,figname);text(0,0.5,'gamma coupling to CSD during theta');text(0,0,'X gamma phase in deg; Y gamma frq in Hz')
    th_spnum=size(gm_cellcwtspect{2,un_indx},2);
    cmax=(th_spnum/20);
    
    for je=2:ch_n-1
        subplot(15,3,49-(je*3))
        set(gca,'FontSize',5)
        if ~isnan(gm_cellcwtspect{je,un_indx})
            for ka=1:size(gm_cellcwtspect{je,un_indx},1)
                gm_phasehist{je,un_indx}(ka,:)=repmat((hist(rad2deg(angle(gm_cellcwtspect{je,un_indx}(ka,:)))+180,[gm_binsize/2:gm_binsize:360-gm_binsize/2])),1,2);
            end; clear ka
            surf([-gm_binsize/2:gm_binsize:720+gm_binsize/2], gm_fraxis, [gm_phasehist{je,un_indx}(:,end) gm_phasehist{je,un_indx} gm_phasehist{je,un_indx}(:,1)],'EdgeColor', 'interp', 'FaceColor', 'interp'); shading flat; colormap jet; caxis([0 2*cmax]);view(2);
            axis([0 720 min(gm_fraxis) max(gm_fraxis)]);
            set(gca, 'YTick',[(min(gm_fraxis)):10:(max(gm_fraxis))]); set(gca, 'XTick',[0:180:720]);
        else gm_phasehist{je}=NaN;
        end
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
    end; clear je cmax
    
    %calculus
    for je=2:ch_n-1
        if ~isnan(gm_cellcwtspect{je,un_indx});
            gm_cellcwtspect_r{un_indx}(je,:)=circ_r(angle(gm_cellcwtspect{je,un_indx}),[],[],2);
            gm_cellcwtspect_m{un_indx}(je,:)=(rad2deg(circ_mean(angle(gm_cellcwtspect{je,un_indx}),[],2)))+180;
            for ka=1:size(gm_cellcwtspect{je,un_indx},1)
                [pe,ze]=circ_rtest(angle(gm_cellcwtspect{je,un_indx}(ka,:)));
                gm_cellcwtspect_p{un_indx}(je,ka)=pe; gm_cellcwtspect_z{un_indx}(je,ka)=ze;
                if pe<=gm_plimit
                    gm_cellcwtspect_sign{un_indx}(je,ka)=1;
                end
                clear pe ze;
            end; clear ka;
        end
    end
    
    %gamma  coupling (during theta); coupling strength spectra
    for je=2:ch_n-1
        subplot(15,3,(50-(je*3)))
        set(gca,'FontSize',5)
        if ~isnan(gm_cellcwtspect{je,un_indx});
            line(gm_cellcwtspect_r{un_indx}(je,:),gm_fraxis);
            axis([0 max(max(gm_cellcwtspect_r{un_indx})*1.1) (min(gm_fraxis)) (max(gm_fraxis)) ]);
            line((gm_cellcwtspect_r{un_indx}(je,:).*gm_cellcwtspect_sign{un_indx}(je,:)), gm_fraxis,'color','red');
            set(gca, 'YTick',[(min(gm_fraxis)):10:(max(gm_fraxis))]);
        end
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
    end; clear je;
    
    %gamma  coupling (during theta); coupling phase spectra
    for je=2:ch_n-1
        subplot(15,3,(51-(je*3)))
        set(gca,'FontSize',5)
        if ~isnan(gm_cellcwtspect{je,un_indx});
            line(gm_cellcwtspect_m{un_indx}(je,:), gm_fraxis);
            line((gm_cellcwtspect_m{un_indx}(je,:))+360, gm_fraxis);
            line((gm_cellcwtspect_m{un_indx}(je,:).*gm_cellcwtspect_sign{un_indx}(je,:)), gm_fraxis,'color','red');
            line((((gm_cellcwtspect_m{un_indx}(je,:))+360).*gm_cellcwtspect_sign{un_indx}(je,:)), gm_fraxis,'color','red');
            axis([0 720 (min(gm_fraxis)) (max(gm_fraxis))]);
            set(gca, 'YTick',[(min(gm_fraxis)):10:(max(gm_fraxis))]); set(gca, 'XTick',[0:90:720]);
        end; clear je;
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
    end
end

%cleanup
clear fl_indx fl_ord ch_def ch_spks fl_name ans figname un_indx       

gmtht.gm_fraxis=gm_fraxis;

gmtht.gm_sptrig_spikenum=th_spnum;
gmtht.gm_sptrig_rawcwt=gm_cellcwtspect;
gmtht.gm_sptrig_phasespect=gm_cellcwtspect_m;
gmtht.gm_sptrig_rspect=gm_cellcwtspect_r;
gmtht.gm_sptrig_zspect=gm_cellcwtspect_z;
gmtht.gm_sptrig_Pspect=gm_cellcwtspect_p;
gmtht.gm_sptrig_signif=gm_cellcwtspect_sign;
gmtht.gm_sptrig_histograms=gm_phasehist;
gmtht.gm_sptrig_plimit=gm_plimit;
gmtht.ch_samplrate=ch_rate;

%saving the relevant data
writename=strcat(bs_name, bs_num, bs_exp, '_a_gma_thtCOUPexpand');
save(writename,'gmtht','-v7.3');




 
    
           
 