%This program generates CSD from LFP array data
clear all
close all

%give the experimenter identifier, the experiment number,m and the experiment type
bs_name='B'; bs_num='186'; bs_unit='d'; bs_typ='sx';

%give the experiment details
dz=0.05;                                            %contact site spacing (in mm)
ch_n=16;                                            %number of channels
%order of channels, separated by space, starting from the ventralmost to the dorsalmost
%ch_ord=[6,11,3,14,1,16,2,15,5,12,4,13,7,10,8,9];                %single shank linear probe
ch_ord=[9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24];       %central shank of the three shank custom probe (pre-ordered)
%ch_ord=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];    %order of channels, separated by space, starting from the ventralmost to the dorsalmost
ch_def=[];                                          %identifier of deffective channels, separated by space; if extreme channels are deffected that decreases the probe site and does not include here
fl_ord=[1];                                         %file extension numbers belonging to the present cell

%give analysis parameters to be used
ch_inv = 1;                                     %channel inversion; write 1 if not needed, write -1 if needed (or a scaling factor if necessary)
sam_rate = 2000;                                    %requested sampling rate in Hz (LFP resampling will be performed if non-matching)
ch_step = 1;                                        %requested step size for the CSD calculus
plotornot = 1;                                	    %should I plot or not the results

%file index stepper
for fl_indx=1:length(fl_ord);
    
    % name of the actually processed file
    fl_name=strcat(bs_name, bs_num, bs_unit, num2str(fl_ord(fl_indx)),'WAW.mat');
    
    %get the length (seconds) of all channels and take the minimum
    fl_length=zeros(size(ch_ord));
    for je=1:ch_n
        ch_name=strcat(bs_name, bs_typ, bs_num, bs_unit,'_Ch', num2str(ch_ord(je))); ch_load=load(num2str(fl_name),num2str(ch_name));
        fl_length(je)=(((ch_load.(ch_name).('length'))-1)*(ch_load.(ch_name).('interval')))+(ch_load.(ch_name).('start'));
    end; clear je;
    fl_length=min(fl_length);
    sam_length=floor(fl_length*sam_rate);
    
    %LFP import channel by channel; resample if necessary
    LFP = zeros(ch_n,sam_length);
    LFP_tbase = [(1/sam_rate):(1/sam_rate):sam_length*(1/sam_rate)];
    for je=1:ch_n
        %load data for the actual channel
        ch_name = strcat(bs_name, bs_typ, bs_num, bs_unit,'_Ch', num2str(ch_ord(je))); ch_load=load(num2str(fl_name),num2str(ch_name));
        LFP_in = ch_inv.*(ch_load.(ch_name).('values'));
        LFP_samint = ch_load.(ch_name).('interval');
        LFP_samst = ch_load.(ch_name).('start');
        
        if isempty(ch_def==ch_ord(je))
            
            %check if the time offset and the sampling interval is the desired
            if LFP_samint == 1/sam_rate && LFP_samst == 1/sam_rate
               %put LFP channels into the final variable
                LFP(je,:) = LFP_in(1:sam_length);
            else
                %resample timeseries and transform it back
                LFP_ts = resample((timeseries(LFP_in, [LFP_samst:LFP_samint:((LFP_samint*(length(LFP_in)-1))+LFP_samst)])),LFP_tbase);
                LFP(je,:) = LFP_ts.Data;
                clear LFP_ts
            end
            clear LFP_samint LFP_samst LFP_in ch_name 
        end
    end;
    clear je ch_load LFP_tbase
    
    %interpolating the deffective channels (linear interpolation from the adjacent channels)
    if ~isempty(ch_def);
        for ka=1:length(ch_def);
            je=find(ch_ord==(ch_def(k)));
            LFP(je,:)=(LFP(je-1,:)+LFP(je+1,:))./2;
        end; clear ka je;
    end
    clear fl_length
    
    %generate the CSD operator matrix
    CSD_matr = zeros(ch_n-(2*ch_step), ch_n);
    for je=1:(ch_n-(2*ch_step))
        CSD_matr(je,je)=1;
        CSD_matr(je,je+ch_step)=-2;
        CSD_matr(je,je+(2*ch_step))=1;
    end; clear je
    
    %CSD calculus
    CSD = (CSD_matr*LFP)./((dz*ch_step)^2);
    CSD_dum(1:ch_step,1:size(CSD,2)) = NaN;
    CSD = [CSD_dum; CSD; CSD_dum];
    
    %writing files
    fl_write=strcat(bs_name, bs_num, bs_unit, num2str(fl_ord(fl_indx)), 'LFP');
    save(fl_write, 'LFP');
    fl_write=strcat(bs_name, bs_num, bs_unit, num2str(fl_ord(fl_indx)), 'CSD');
    save(fl_write, 'CSD');
    
    if plotornot==1
        for je=1:ch_n
            LFPplot(je,:)=LFP(je,1:5000)+(0.7*je);
        end; clear je;
        for je=1:(ch_n-2)
            CSDplot(je,:)=CSD(je,1:5000)+(100*je);
        end; clear je;
        figure
        plot(LFPplot');xlim([1 2000]);
        figure
        plot(CSDplot');xlim([1 2000]);
    end
    clear CSD LFP CSD_dum CSDplot LFPplot CSD_matr fl_length fl_name fl_write LFP_tbase sam_length
end

   

    
            