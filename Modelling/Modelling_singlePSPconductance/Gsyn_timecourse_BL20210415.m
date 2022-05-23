%exploring the experimental reationship between 10-90% risetime and tau
%rise in the equation shows that tau of 4.55 corresponds to 10-90%
%rise-time of 4.2 ms (a tau decay of 20 ms was used)

tdec=35
tris=4.24
%for je=1:300
%    tris=3.00+(je/100);
    
    for stp=1:30001
        tm=(stp-1)/100;
        tmscale(stp)=tm;
        gsyn(stp)=(1-(exp((-1*tm)/tris)))*(exp((-1*tm)/tdec));
    end
    [pkval,pktim]=max(gsyn);
    gsynE=gsyn(1:max(pktim));
    [~,inds]=find(gsynE<(0.1*pkval));
    ind10=max(inds);
    [~,inds]=find(gsynE>(0.9*pkval));
    ind90=min(inds);
    rtm=(ind90-ind10)/100;
    
%     iksz(je)=tris; ipsz(je)=rtm; pktm(je)=pktim/100;
% end