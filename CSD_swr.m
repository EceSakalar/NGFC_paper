function [CSD_rips] = CSD_swr(CSD, ch, swr, rp)
%%CSD%%
%preallocate variable
for je=1:ch_n
    CSD_rips{je}=zeros(length(rp.fraxis),sum(n_tbins));
end
CSD_ripsN=CSD_rips;
clear je
CSD_segavg=zeros(ch.n,sum(n_tbins));

for je=(ch.step+1):(ch.n-ch.step)  %CSD channel stepper
    if isempty(ch.def(ch.def==(je))) %exclude defected channels
        
        CSDact=CSD(je,:)*-1;
        clear CSD csd_name;
        
        if isnan(CSDact)
            CSD_rips{je}=NaN;
            CSD_segavg(je,:)=NaN;
        else
            for ka=1:size(swr.pers,1)
                sta_rip=round(ch_rate*swr.pers(ka,1));
                end_rip=round(ch_rate*swr.pers(ka,2));
                bin_siz=(end_rip-sta_rip)/n_tbins(2);
                sta_wlt=floor(sta_rip-(bin_siz*(n_tbins(1)*1.4)));
                end_wlt=ceil(end_rip+(bin_siz*(n_tbins(3)*1.4)));
                csd_seg=CSDact(sta_wlt:end_wlt);
                wlt_rip=conj(cwt(csd_seg, rp_scsta:rp.scint:rp.scfin, rp.wavelet{1}));
                sta_rip=sta_rip-sta_wlt+1;
                end_rip=end_rip-sta_wlt+1;
                bin_1st=sta_rip-(bin_siz*n_tbins(1));
                amp_rip=zeros(length(rp_fraxis),sum(n_tbins));
                for el=0:sum(n_tbins)-1
                    wlt_el=abs(wlt_rip(:,(ceil(bin_1st+(el*bin_siz))):(floor(bin_1st+((el+1)*bin_siz)))));
                    csd_el=csd_seg((ceil(bin_1st+(el*bin_siz))):(floor(bin_1st+((el+1)*bin_siz))));
                    amp_rip(:,el+1)=mean(wlt_el,2); clear wlt_el
                    amp_csd(el+1)=mean(csd_el); clear csd_el
                end
                CSD_rips{je}=CSD_rips{je}+amp_rip; clear amp_rip
                CSD_segavg(je,:)=CSD_segavg(je,:)+amp_csd;
                clear amp_csd
            end
            clear sta_rip end_rip bin_siz sta_wlt end_wlt wlt_rip bin_1st ka csd_seg CSDact
        end
    else
        CSD_rips{je}=NaN;
        CSD_ripsN{je}=NaN;
    end
end
clear je swr.pers

end

