function [] = thUntCoupl_Plt(bs, thcell)
%PLOTTING STATISTICS
%plotting figures
figindex=0;

for un_indx=1:length(un_chans)
    figindex_T=ceil(un_indx/6);
    if figindex_T~=figindex
        figindex=figindex+1;
        figname=strcat(bs.name, bs.num, bs.exp, num2str(figindex),'_thetamodulation');
        figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
        plotindex=1;
    else
        plotindex=((un_indx-((figindex-1)*6))*3)-2;
    end
    if isempty(thcell.coupvect_lfp{un_indx})
        subplot(6,3,plotindex);
        axis off;
        text(0, 1, 'no data for this cell');
    else
        %textual results
        subplot (6,3,plotindex);
        axis off;
        actualtext=strcat(bs.name, bs.num, bs.exp, 'sh', num2str(un_shnk(un_indx)),'u',num2str(un_nums(un_indx))); 
        text(0, 1, actualtext);
        %number of theta cycles
        actualtext=strcat('Number of theta cycles analysed: ', num2str(thcell.cyccount(un_indx))); 
        text(0, 0.9, actualtext);
        %number of active theta cycles
        th_acyccount=length(thcell.spikespercycles{un_indx});
        th_acycperc=(th_acyccount/thcell.cyccount(un_indx))*100; 
        actualtext=strcat('of these active: ',num2str(th_acyccount),' (',num2str(th_acycperc),'%)'); 
        text(0, 0.8, actualtext);
        %%mean number of APs in active theta cycles
        actualtext=strcat('# of spikes/cycle (SD):',num2str(mean(thcell.spikespercycles{un_indx})),' (',num2str(std(th_cycspkcount{un_indx})),')'); 
        text(0, 0.7, actualtext);
        %mean AP frequency
        th_spkrate=length(thcell.coupvect_lfp{un_indx})/th_time(un_indx); 
        actualtext=strcat('spike frequency under theta:', num2str(th_spkrate), ' s-1'); 
        text(0, 0.6, actualtext);
        %mean vector legth
        actualtext=strcat('mean vector length:', num2str(thcell.coupvect_lfp_r(un_indx))); 
        text(0, 0.5, actualtext);
        %mean phase angle
        actualtext=strcat('mean vector angle:', num2str(thcell.coupvect_lfp_m(un_indx))); 
        text(0, 0.4, actualtext);
        % angular deviation and circular SD
        s=(180/pi)*(sqrt(2*(1-thcell.coupvect_lfp_r(un_indx)))); 
        so=(180/pi)*(sqrt(-2*(log(thcell.coupvect_lfp_r(un_indx)))));
        actualtext=strcat('Angular deviation (s): ', num2str(s)); 
        text(0, 0.3, actualtext);
        actualtext=strcat('Circular standard deviation (so): ', num2str(so)); 
        text(0, 0.2, actualtext);
        %Rayleigh statistics
        n=length(thcell.coupvect_lfp{un_indx});
        R=n*thcell.coupvect_lfp_r(un_indx); 
        z=R.^2/n; 
        P=exp((sqrt(1+(4*n)+(4*(n.^2-R.^2))))-(1+(2*n)));
        text(0, 0, 'Rayleigh test for the uniformity of distribution');
        actualtext=strcat ('z=', num2str(z),' P=', num2str(P),'n=', num2str(n)); 
        text(0, 0.1, actualtext);
        %expand the graph
        pozi=get(gca,'Position'); 
        pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; 
        set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
        
        subplot (6,3,plotindex+1);
        cycbintimes=[repmat(sum(cycdur{un_indx}(:,1),1)./10, 1, th_nbins/2) repmat(sum(cycdur{un_indx}(:,2),1)./10, 1, th_nbins/2)];
        line([binscale binscale+360], repmat(thcell.coupvect_lfp_h(un_indx,:)./(cycbintimes),1,2));
        x=(0:1:720); 
        y=(cosd(x+180)+1)*(max(thcell.coupvect_lfp_h(un_indx,:)./max(cycbintimes))).*0.5; 
        line(x,y, 'Color', 'red');
        axis([0 720 0 max(thcell.coupvect_lfp_h(un_indx,:)./(th_time(un_indx)/20))*1.2]);
        set(gca, 'Xtick', 0:90:720); 
        xlabel('theta phase (deg)'); 
        ylabel('spike rate (s-1)')
        %expand the graph
        pozi=get(gca,'Position'); 
        pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; 
        set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
        
        %phase plot for all spikes
        subplot (6,3,plotindex+2);
        plot(thcell.coupvect_lfp{un_indx}, (1:1:length(thcell.coupvect_lfp{un_indx})), 'b. ');
        axis([0 360 0 length(thcell.coupvect_lfp{un_indx})]); 
        xlabel('phase (deg)'); 
        ylabel('spike #');
        %expand the graph
        pozi=get(gca,'Position'); 
        pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; 
        set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
        
    end
end
clear th_nbins;
end

