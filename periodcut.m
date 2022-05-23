%this functions removes all (3rd argument 0) or preserves only (3rd
%argument 1) the elements of vector (1st argument) if they fall within 
%the limits defined by any of the lines of period matrix (2nd argument)

function [cutevents]=periodcut(origevents,periods,argum)

cutevents=[];

if argum(1)==0
    if isempty(periods)==0
        for k=1:size(periods, 1);
            origevents(find(origevents>(periods(k,1))&origevents<(periods(k,2))))=[];
        end
        cutevents=origevents;
    else cutevents=origevents;
    end
end

if argum(1)==1
    if isempty(periods)==0
        for k=1:size(periods, 1);
            cutevents_temp=origevents(find(origevents>=(periods(k,1))&origevents<(periods(k,2))));
            cutevents=[cutevents; cutevents_temp]; clear cutevents_temp;
        end
    else cutevents=[];
    end
end
    

    