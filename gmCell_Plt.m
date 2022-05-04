function [ output_args ] = gmCell_Plt(gm)
%plotting results gamma coupling (during theta)
%calculus
for un_indx=1:length(un_chans)
    gm_cellcwtspect_r{un_indx}(1:ch_n,1:length(gm.fraxis))=NaN;
    gm_cellcwtspect_m{un_indx}(1:ch_n,1:length(gm.fraxis))=NaN;
    gm_cellcwtspect_p{un_indx}(1:ch_n,1:length(gm.fraxis))=NaN;
    gm_cellcwtspect_z{un_indx}(1:ch_n,1:length(gm.fraxis))=NaN;
    gm_cellcwtspect_sign{un_indx}(1:ch_n,1:length(gm.fraxis))=NaN;
end

for un_indx=1:length(un_chans)
    %histograms
    figname=strcat(bs.name, bs.num, bs.exp, ' shk', num2str(un_shank(un_indx)),'u',num2str(un_nums(un_indx)));
    figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
    subplot(15,3,1);
    axis off;
    text(0,1,figname);
    text(0,0.5,'gamma coupling to CSD during theta');
    text(0,0,'X gamma phase in deg; Y gamma frq in Hz')
    th_spnum=size(gm_cellcwtspect{2,un_indx},2);
    cmax=(th_spnum/20);
    
    for je=2:ch_n-1
        subplot(15,3,49-(je*3))
        set(gca,'FontSize',5)
        if ~isnan(gm_cellcwtspect{je,un_indx})
            for ka=1:size(gm_cellcwtspect{je,un_indx},1)
                gm_phasehist{je,un_indx}(ka,:)=repmat((hist(rad2deg(angle(gm_cellcwtspect{je,un_indx}(ka,:)))+180,[gm_binsize/2:gm_binsize:360-gm_binsize/2])),1,2);
            end; clear ka
            surf([-gm_binsize/2:gm_binsize:720+gm_binsize/2], gm.fraxis, [gm_phasehist{je,un_indx}(:,end) gm_phasehist{je,un_indx} gm_phasehist{je,un_indx}(:,1)],'EdgeColor', 'interp', 'FaceColor', 'interp'); 
            shading flat; 
            colormap jet; 
            caxis([0.5*cmax 1.5*cmax]);
            view(2);
            axis([0 720 min(gm.fraxis) max(gm.fraxis)]);
            set(gca, 'YTick',[(min(gm.fraxis)):10:(max(gm.fraxis))]); 
            set(gca, 'XTick',[0:180:720]);
        else
            gm_phasehist{je}=NaN;
        end
        %expand the graph
        pozi=get(gca,'Position'); 
        pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; 
        set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
    end
    clear je cmax
    
    %calculus
    for je=2:ch_n-1
        if ~isnan(gm_cellcwtspect{je,un_indx})
            for ka=1:size(gm_cellcwtspect{je,un_indx},1)
                gm_cellcwtspect_m{un_indx}(je,ka)=(rad2deg(circ_mean(angle(gm_cellcwtspect{je,un_indx}(ka,:)))))+180;
                gm_cellcwtspect_r{un_indx}(je,ka)=circ_r(angle(gm_cellcwtspect{je,un_indx}(ka,:)));
                [pe,ze]=circ_rtest(angle(gm_cellcwtspect{je,un_indx}(ka,:)));
                gm_cellcwtspect_p{un_indx}(je,ka)=pe; 
                gm_cellcwtspect_z{un_indx}(je,ka)=ze;
                if pe<=gm_plimit
                    gm_cellcwtspect_sign{un_indx}(je,ka)=1;
                end
                clear pe ze;
            end 
            clear ka;
        end
    end
    
    %gamma  coupling (during theta); coupling strength spectra
    for je=2:ch_n-1
        subplot(15,3,(50-(je*3)))
        set(gca,'FontSize',5)
        if ~isnan(gm_cellcwtspect{je,un_indx})
            line(gm_cellcwtspect_r{un_indx}(je,:),gm.fraxis);
            axis([0 max(max(gm_cellcwtspect_r{un_indx})*1.1) (min(gm.fraxis)) (max(gm.fraxis)) ]);
            line((gm_cellcwtspect_r{un_indx}(je,:).*gm_cellcwtspect_sign{un_indx}(je,:)), gm.fraxis,'color','red');
            set(gca, 'YTick',[(min(gm.fraxis)):10:(max(gm.fraxis))]);
        end
        %expand the graph
        pozi=get(gca,'Position'); 
        pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; 
        set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
    end 
    clear je;
    
    %gamma  coupling (during theta); coupling phase spectra
    for je=2:ch_n-1
        subplot(15,3,(51-(je*3)))
        set(gca,'FontSize',5)
        if ~isnan(gm_cellcwtspect{je,un_indx})
            line(gm_cellcwtspect_m{un_indx}(je,:), gm.fraxis);
            line((gm_cellcwtspect_m{un_indx}(je,:))+360, gm.fraxis);
            line((gm_cellcwtspect_m{un_indx}(je,:).*gm_cellcwtspect_sign{un_indx}(je,:)), gm.fraxis,'color','red');
            line((((gm_cellcwtspect_m{un_indx}(je,:))+360).*gm_cellcwtspect_sign{un_indx}(je,:)), gm.fraxis,'color','red');
            axis([0 720 (min(gm.fraxis)) (max(gm.fraxis))]);
            set(gca, 'YTick',[(min(gm.fraxis)):10:(max(gm.fraxis))]); 
            set(gca, 'XTick',[0:90:720]);
        end 
        clear je;
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
    end
    
    %saving figure
    fg_name=strcat('fig',bs.name, bs.num, bs.exp,'GMAcellTHT_shk', num2str(un_shank(un_indx)),'u',num2str(un_nums(un_indx)));
    saveas(gcf,fg_name)
    clear fg_name
end

end

