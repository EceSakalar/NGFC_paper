function [CSD] = CSD_cal(LFP, ch)
%CSD calculation from LFP
%generate the CSD operator matrix

CSD_matr = zeros(ch.n-(2*ch.step), ch.n);

for je=1:(ch.n-(2*ch.step))
    CSD_matr(je,je)=1;
    CSD_matr(je,je+ch.step)=-2;
    CSD_matr(je,je+(2*ch.step))=1;
end
clear je

%CSD calculus
CSD = (CSD_matr*LFP)./((ch.dz*ch.step)^2);
CSD_dum(1:ch.step,1:size(CSD,2)) = NaN;
CSD = [CSD_dum; CSD; CSD_dum];

end

