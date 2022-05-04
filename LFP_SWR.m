function [LFP_rips] = LFP_SWR(LFP, ch, swr, rp)
%This program analyses the SWR coupling and ripple oscillatory coupling relative to LFP and CSD traces

%SWR counter reset
SWR_no=0;
%%LFP%%
%preallocate variables
for je=1:ch.n
    LFP_rips{je}=zeros(length(rp.fraxis),sum(n_tbins));
end
LFP_ripsN=LFP_rips;
clear je
LFP_segavg=zeros(ch.n,sum(n_tbins));

for je=1:ch.n  %LFP channel stepper
    if isempty(ch.def(ch.def==(je))) %exclude defected channels
        
        %reading LFP for ripple LFP characterisation
        LFPact=LFP(je,:);
        for ka=1:size(swr.pers,1)
            if je==1
                SWR_no=SWR_no+1;
            end
            sta_rip=round(ch_rate*swr.pers(ka,1));
            end_rip=round(ch_rate*swr.pers(ka,2));
            bin_siz=(end_rip-sta_rip)/n_tbins(2);
            sta_wlt=floor(sta_rip-(bin_siz*(n_tbins(1)*1.4)));
            end_wlt=ceil(end_rip+(bin_siz*(n_tbins(3)*1.4)));
            lfp_seg=LFPact(sta_wlt:end_wlt);
            wlt_rip=conj(cwt(lfp_seg,rp_scsta:rp_scint:rp_scfin, rp_wavelet{1}));
            sta_rip=sta_rip-sta_wlt+1;
            end_rip=end_rip-sta_wlt+1;
            bin_1st=sta_rip-(bin_siz*n_tbins(1));
            amp_rip=zeros(length(rp_fraxis),sum(n_tbins));
            for el=0:sum(n_tbins)-1
                wlt_el=abs(wlt_rip(:,(ceil(bin_1st+(el*bin_siz))):(floor(bin_1st+((el+1)*bin_siz)))));
                lfp_el=lfp_seg((ceil(bin_1st+(el*bin_siz))):(floor(bin_1st+((el+1)*bin_siz))));
                amp_rip(:,el+1)=mean(wlt_el,2); clear wlt_el
                amp_lfp(el+1)=mean(lfp_el); clear lfp_el
            end
            clear el
            LFP_rips{je}=LFP_rips{je}+amp_rip; clear amp_rip
            LFP_segavg(je,:)=LFP_segavg(je,:)+amp_lfp; clear amp_lfp
        end
        clear sta_rip end_rip bin_siz sta_wlt end_wlt wlt_rip bin_1st ka lfp_seg LFPact
        clear swr.pers
    else
        LFP_rips{je}=NaN;
        LFP_ripsN{je}=NaN;
    end
end
clear je

end

