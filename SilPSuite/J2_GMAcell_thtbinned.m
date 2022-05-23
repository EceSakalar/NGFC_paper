%This program analyses the gamma oscillatory coupling relative to CSD
%traces (spectral aproach; multichannel LFP and CSD traces), binned for
%different theta phases to offset the co-modulation of unit spikes and
%gamma oscillations by theta cycles

%the data are read from results files created earlier by the theta and
%gamma scripts

clear all; close all
bs_name='B'; bs_num='172'; bs_exp='c';                                     %unit basic data entry

%give the number of theta bins to be used and the min number of spikes in
%bins to allow processing
tht_binnum=6; tht_binmin=200; tht_binsize=360/tht_binnum;

%give gamma statistic data
gma_plimit=0.05;

%load the theta phases for all spikes
flT_name=strcat(bs_name, bs_num, bs_exp, '_a_tht_cellcoupling.mat'); load (flT_name); clear flT_name
tht_phas=thcell.coupvect_lfp;

%unit properties
un_nums=[33 40 51 53 54 65 67 10 19];                               %unit numbers (in kwik file)
un_shank =[repmat([1],1,7) repmat([2],1,2)];                        %shanks from which unit was isolated

%assign spikes to theta bins
for je=1:size(tht_phas,2);
    tht_bins{1,je}=(round(tht_phas{1,je}./tht_binsize))+1;
    tht_bins{1,je}(tht_bins{1,je}==(tht_binnum+1))=1;
    bin_spkn{1,je}=min(hist(tht_bins{1,je},1:tht_binnum));
end; clear je

%load the gamma instantaneous wavelets for all the spikes
flG_name=strcat(bs_name, bs_num, bs_exp, '_a_gma_thtCOUP'); load(flG_name); clear flG_name
gma_icwt=gmtht.gm_sptrig_rawcwt;
gma_fraxis=gmtht.gm_fraxis;

%take gamma phase samples with biased sample numbers
for je=1:size(tht_phas,2);
    if bin_spkn{je}<tht_binmin
        for ka=1:size(gma_icwt,1)
            for el=1:tht_binnum;
                gma_phas{ka,je,el}=[];
            end; clear el
        end; clear ka
    else
        for el=1:tht_binnum;
            %create a sampler for the theta bin with bin_spkn random spikes
            [~, bin_sampler]=sort(rand(size(tht_bins{je}(tht_bins{je}==el),2),1));
            for ka=1:size(gma_icwt,1)
                if isempty(gma_icwt{ka,je});
                    gma_phas{ka,je,el}=[];
                else
                    gma_bcwt{ka,je,el}=gma_icwt{ka,je}(:,tht_bins{je}==el);
                    gma_phas{ka,je,el}=(angle(gma_bcwt{ka,je,el}(:,sort(bin_sampler(1:bin_spkn{je})))))+pi;
                end
            end; clear ka
        end; clear el bin_sampler
    end
end; clear je

%------------------------------generating plots----------------------------
%take gamma phase samples with biased numbers and generate stat data
for je=1:size(tht_phas,2);
    for ka=1:size(gma_icwt,1);
        if isempty(gma_phas{ka,je})
            gma_phRb{ka,je}=[];
            gma_phMb{ka,je}=[];
            gma_phPb{ka,je}=[];
            gma_phSb{ka,je}=[];
            for el=1:tht_binnum;
                gma_phHb{ka,je,el}=[];
            end
        else
            for el=1:tht_binnum
                gma_phRb{ka,je}(:,el)=circ_r(gma_phas{ka,je,el}, [],[],2);
                gma_phMb{ka,je}(:,el)=rad2deg(circ_mean(gma_phas{ka,je,el}')); gma_phMb{ka,je}(gma_phMb{ka,je}<0)=gma_phMb{ka,je}(gma_phMb{ka,je}<0)+360;
                for em=1:size(gma_phas{ka,je,el},1)
                    [pe,ze]=circ_rtest(gma_phas{ka,je,el}(em,:));
                    gma_phPb{ka,je}(em,el)=pe;
                    if pe<=gma_plimit
                        gma_phSb{ka,je}(em,el)=1;
                    else
                        gma_phSb{ka,je}(em,el)=NaN;
                    end; clear pe ze;
                end; clear em
            end; clear el;
        end
    end; clear ka
end; clear je
                
%-------------------------------plotting--------------------------------------
for je=1:size(tht_phas,2);
    if ~isempty(gma_phRb{2,je})
        figname=strcat(bs_name, bs_num, bs_exp, ' shk', num2str(un_shank(je)),'u',num2str(un_nums(je)));
        figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
        subplot(15,2,1);
        text(0,1,figname);text(0,0.5,'gamma coupling to CSD during theta bins');text(0,0,'X gamma phase in deg; Y gamma frq in Hz')
        axis off
        
        %theta bins
        subplot(15,2,2);
        colormap jet
        AlCols=colormap;
        LiCols=AlCols((round(((1:tht_binnum)./tht_binnum).*64)),:); clear AlCols;
        for el=1:tht_binnum
            tht_hist(1,el)=(el-1)*tht_binsize;
            tht_hist(2,el)=length(tht_bins{je}(tht_bins{je}==el));
        end; clear el;
        tht_histx=[tht_hist(1,:) tht_hist(1,:)+360 tht_hist(1,:)+720];
        tht_histy=[tht_hist(2,:) tht_hist(2,:) tht_hist(2,:)];
        HiCols=[LiCols; LiCols; LiCols];
        scatter(tht_histx, tht_histy, 20, HiCols, 'filled')
        axis([0 720 0 max(tht_histy)])
        set(gca, 'XTick',[0:180:720])
            
        %calculate max R
        for ka=2:15
            Rmaxs(ka)=max(max(gma_phRb{ka,je}));
        end
        
        for ka=2:15
            %mean vector length spectra
            subplot(15,2,31-((ka-1)*2))
            set(gca,'FontSize',5)
            for el=1:tht_binnum
                line(gma_phRb{ka,je}(:,el),gma_fraxis, 'Color', [0.8 0.8 0.8]);
                if je==1; hold on; end;
                line(((gma_phRb{ka,je}(:,el)).*(gma_phSb{ka,je}(:,el))),gma_fraxis, 'Color', LiCols(el,:));
                if je==tht_binnum; hold off; end;
            end; clear el
            axis([0 max(Rmaxs) (min(gma_fraxis)) (max(gma_fraxis))]);
            %expand the graph
            pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi)
            set(gca, 'Clipping', 'off')
            
            %mean phase spectra
            subplot(15,2,32-((ka-1)*2))
            set(gca,'FontSize',5)
            for el=1:tht_binnum
                line(gma_phMb{ka,je}(:,el),gma_fraxis, 'Color', [0.8 0.8 0.8]);
                if je==1; hold on; end;
                line(gma_phMb{ka,je}(:,el)+360,gma_fraxis, 'Color', [0.8 0.8 0.8]);
                line((gma_phMb{ka,je}(:,el).*(gma_phSb{ka,je}(:,el))),gma_fraxis', 'Color', LiCols(el,:));
                line(((gma_phMb{ka,je}(:,el)).*(gma_phSb{ka,je}(:,el)))+360, gma_fraxis', 'Color', LiCols(el,:));
                if je==tht_binnum; hold off; end;
            end;clear el
            axis([0 720 (min(gma_fraxis)) (max(gma_fraxis))]);
            %expand the graph
            pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi)
            set(gca, 'Clipping', 'off')
        end; clear ka
    end
end; clear je


%cleanup
clear gmtht HiCols pozi Rmaxs thcell figname un_nums un_shank tht_histx tht_histy tht_hist; 

gmthtbin.gm_fraxis=gma_fraxis;
gmthtbin.gm_rawcwt=gma_icwt;
gmthtbin.gm_phspect_select=gma_phas;
gmthtbin.gm_rspect=gma_phRb;
gmthtbin.gm_mspect=gma_phMb;
gmthtbin.gm_pspect=gma_phPb;
gmthtbin.gm_sspect=gma_phSb;
gmthtbin.gm_plimit=gma_plimit;

gmthtbin.th_bins=tht_bins;
gmthtbin.th_phases=tht_phas;
gmthtbin.th_binnum=tht_binnum;
gmthtbin.th_minspknums=bin_spkn;

%saving the relevant data
writename=strcat(bs_name, bs_num, bs_exp, '_a_gma_thtbinnedCOUP');
save(writename,'gmthtbin','-v7.3');




 
    
           
 