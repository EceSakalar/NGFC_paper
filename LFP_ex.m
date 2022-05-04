function [LFP] = LFP_ex(fl_name, ch)
%LFP extraction from .mat files exported from Spike2 recordings

%get the length (seconds) of all channels and take the minimum
fl_length=zeros(size(ch.ord));
for je=1:ch.n
    ch_name=strcat('Lin1_', num2str(ch.ord(je).','%02d'));
    ch_load=load(num2str(fl_name),num2str(ch_name));
    fl_length(je)=(((ch_load.(ch_name).('length'))-1)*(ch_load.(ch_name).('interval')))+(ch_load.(ch_name).('start'));
end; clear je;
fl_length=min(fl_length);
sam_length=floor(fl_length*ch.SampRate);

%LFP import channel by channel; resample if necessary
LFP = zeros(ch.n,sam_length);
LFP_tbase = [(1/ch.SampRate):(1/ch.SampRate):sam_length*(1/ch.SampRate)];
for je=1:ch.n
    %load data for the actual channel
    ch_name=strcat('Lin1_', num2str(ch.ord(je).','%02d'));
    ch_load=load(num2str(fl_name),num2str(ch_name));
    LFP_in = ch.inv.*(ch_load.(ch_name).('values'));
    LFP_samint = ch_load.(ch_name).('interval');
    LFP_samst = ch_load.(ch_name).('start');
    
    if isempty(ch.def==ch.ord(je))
        
        %check if the time offset and the sampling interval is the desired
        if LFP_samint == 1/ch.SampRate && LFP_samst == 1/ch.SampRate
            %put LFP channels into the final variable
            LFP(je,:) = LFP_in(1:sam_length);
        else
            %resample timeseries and transform it back
            LFP_ts = resample((timeseries(LFP_in, [LFP_samst:LFP_samint:((LFP_samint*(length(LFP_in)-1))+LFP_samst)])),LFP_tbase);
            LFP(je,:) = LFP_ts.Data;
            clear LFP_ts
        end
        
        %check the amplitude unit, end value should be in mV
        if ch_load.(ch_name).('units')=='uV'
            LFP(je,:)=LFP(je,:)/1000;
        else
            LFP(je,:)=LFP(je,:);
        end
        
        clear LFP_samint LFP_samst LFP_in ch_name
    end
    clear je ch_load 
    
    %interpolating the deffective channels (linear interpolation from the adjacent channels)
    if ~isempty(ch.def)
        for ka=1:length(ch.def)
            je=find(ch.ord==(ch.def(ka)));
            LFP(je,:)=(LFP(je-1,:)+LFP(je+1,:))./2;
        end 
        clear ka je;
    end
    clear fl_length
end
clear LFP_tbase

end

