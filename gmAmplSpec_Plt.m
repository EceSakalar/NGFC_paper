function [] = gmAmplSpec_Plt(bs, ch_n, gmtht)
%gamma amplitude spectra
figname=strcat(bs.name, bs.num, bs.unit,' gamma amplitude during theta');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
subplot(4,4,1);
axis([0 1 0 4]); 
axis off;
text(0,3,strcat(bs.name,bs.num,bs.unit));
text(0,2,'gamma CSD ampl');
text(0,1,'during theta')
for je=2:ch_n-2 
    smax_t(je)=max(gmtht.gm_osc_amplspectMEAN{je}); 
end 

smax=max(smax_t);
clear smax_t;
 
for je=2:ch_n-1
    subplot(4,4,17-je)
    if ~isnan(gmtht.gm_osc_amplspectMEAN{je})
        xlabel('CSD amplitude (mv/mm2)'); 
        ylabel ('frequency (Hz)');
        set(gca, 'YTick',[(min(gmtht.gm_fraxis)):10:(max(gmtht.gm_fraxis))]);
        plot(gmtht.gm_osc_amplspectMEAN{je}, gmtht.gm_fraxis, gmtht.gm_osc_amplspectMEDIAN{je}, gmtht.gm_fraxis, 'r', gmtht.gm_osc_ampldtdvs{je}, gmtht.gm_fraxis, '--b')
        axis([0 smax*1.1 (min(gmtht.gm_fraxis)) (max(gmtht.gm_fraxis))]);
    end
end
clear je cmax

end