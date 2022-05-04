function [] = gmThMod_Plt(bs, ch_n, gmtht)
%theta-gamma coordination at the reference electrode
figname=strcat(bs.name, bs.num, bs.unit,' gamma modulation by theta');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
subplot(4,4,1);
axis([0 1 0 4]); 
axis off;
text(0,3,strcat(bs.name,bs.num,bs.unit));
text(0,2,'gamma CSD');
text(0,1,'by theta')
for je=1:ch_n
    if isnan(gmtht.gm_osc_thtmod{je})
        cmax_t(je)=0;
    else
        cmax_t(je)=max(max(abs(gmtht.gm_osc_thtmod{je})));
    end
end 
cmax=max(cmax_t); 
clear cmax_t;

 
for je=2:ch_n-1
    subplot(4,4,17-je)
    if ~isnan(gmtht.gm_osc_thtmod{je})
        surf([-th_binsize/2:th_binsize:360+th_binsize/2], gmtht.gm_fraxis, [gmtht.gm_osc_thtmod{je}(:,end) gmtht.gm_osc_thtmod{je} gmtht.gm_osc_thtmod{je}(:,1)], 'EdgeColor', 'interp', 'FaceColor', 'interp'); 
        shading flat; 
        axis tight; 
        colormap jet; 
        caxis([-1*cmax cmax]); 
        view(2);
        axis([0 360 min(gmtht.gm_fraxis) max(gmtht.gm_fraxis)]);
        xlabel('theta phase (deg)'); 
        ylabel ('frequency (Hz)'); 
        set(gca, 'YTick',[(min(gmtht.gm_fraxis)):10:(max(gmtht.gm_fraxis))]); 
        set(gca, 'XTick',[0:180:360]);
    end
end; 
clear je cmax

end
