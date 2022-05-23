%This program generates CSD from LFP array data for Multishank experiments
clear all
close all

%give the experimenter identifier, the experiment number,m and the experiment type
bs_name='IV'; bs_num='004'; bs_unit='a'; bs_typ='sp';

%give the experiment details
dz=0.02;                                            %contact site spacing (in mm)
ch_n=[15 15 15 15 15 15 15 15 3];                     %number of channels for each shank
ch_ord={[1:15], [16:30], [31:45], [46:60], [61:75], [76:90], [91:105], [106:120], [123:125]};         %order of channels, separated by space within shanks, semicolon between shanks starting from the ventralmost to the dorsalmost
ch_def={[],[],[],[],[3],[],[14],[11],[]};                  %identifier of deffective channels, separated by space; if extreme channels are deffected that decreases the probe site and does not include here
fl_ord=[1];                                         %file extension numbers belonging to the present cell

%give analysis parameters to be used
ch_inv = 1;                                         %channel inversion; write 1 if not needed, write -1 if needed (or a scaling factor if necessary)
sam_rate = 10000;                                   %requested sampling rate in Hz (LFP resampling will be performed if non-matching)
ch_step = 1;                                        %requested step size for the CSD calculus
plotornot = 1;                                	    %should I plot or not the results

%shank stepper
for sh_indx=1:length(ch_n)
    
    %file index stepper
    for fl_indx=1:length(fl_ord)
        
        % name of the actually processed file
        fl_name=strcat(bs_name, bs_num, bs_unit, num2str(fl_ord(fl_indx)),'WAW.mat');
        
        %get the length (seconds) of all channels and take the minimum
        fl_length{sh_indx}=zeros(size(ch_ord{sh_indx}));
        for je=1:ch_n(sh_indx)
            ch_name=strcat(bs_name, bs_typ, bs_num, bs_unit,'_Ch', num2str(ch_ord{sh_indx}(je))); ch_load=load(num2str(fl_name),num2str(ch_name));
            fl_inval=ch_load.(ch_name).('interval');
            fl_start=ch_load.(ch_name).('start');
            if fl_start==0; fl_start=fl_inval;end
            fl_length{sh_indx}(je)=(((ch_load.(ch_name).('length'))-1)*(fl_inval))+(fl_start);
        end; clear je fl_inval fl_start ;
        fl_length=min(fl_length{sh_indx});
        sam_length=floor(fl_length*sam_rate);
        
        %LFP import channel by channel; resample if necessary
        LFP = zeros(ch_n(sh_indx),sam_length);
        LFP_tbase = [(1/sam_rate):(1/sam_rate):sam_length*(1/sam_rate)];
        for je=1:ch_n(sh_indx)
            %load data for the actual channel
            ch_name = strcat(bs_name, bs_typ, bs_num, bs_unit,'_Ch', num2str(ch_ord{sh_indx}(je))); ch_load=load(num2str(fl_name),num2str(ch_name));
            LFP_in = ch_inv.*(ch_load.(ch_name).('values'));
            LFP_samint = ch_load.(ch_name).('interval');
            LFP_samst = ch_load.(ch_name).('start');
            if LFP_samst==0; LFP_samst=LFP_samint;end
            
            if isempty(find(ch_def{sh_indx}==ch_ord{sh_indx}(je), 1))
                
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
        end
        clear je ch_load LFP_tbase sam_length
        
        %interpolating the deffective channels (linear interpolation from the adjacent channels)
        if ~isempty(ch_def{sh_indx})
            for ka=1:length(ch_def{sh_indx})
                je=ch_def{sh_indx}(ka);
                LFP(je,:)=(LFP(je-1,:)+LFP(je+1,:))./2;
            end; clear ka je;
        end
        clear fl_length
        
        %generate the CSD operator matrix
        CSD_matr = zeros(ch_n(sh_indx)-(2*ch_step), ch_n(sh_indx));
        for je=1:(ch_n(sh_indx)-(2*ch_step))
            CSD_matr(je,je)=1;
            CSD_matr(je,je+ch_step)=-2;
            CSD_matr(je,je+(2*ch_step))=1;
        end; clear je
        
        %CSD calculus
        CSD = (CSD_matr*LFP)./((dz*ch_step)^2);
        CSD_dum(1:ch_step,1:size(CSD,2)) = NaN;
        CSD = [CSD_dum; CSD; CSD_dum];
        
        %writing files
        fl_write=strcat(bs_name, bs_num, bs_unit, num2str(fl_ord(fl_indx)),'sh',num2str(sh_indx),'LFP');
        save(fl_write, 'LFP','-v7.3');
        fl_write=strcat(bs_name, bs_num, bs_unit, num2str(fl_ord(fl_indx)),'sh',num2str(sh_indx),'CSD');
        save(fl_write, 'CSD','-v7.3');
        
        if plotornot==1
            for je=1:ch_n(sh_indx)
                LFPplot(je,:)=LFP(je,1:5000)+(0.7*je);
            end; clear je;
            for je=1:(ch_n(sh_indx)-2)
                CSDplot(je,:)=CSD(je,1:5000)+(100*je);
            end; clear je;
            figure
            plot(LFPplot');xlim([1 2000]);
            figure
            plot(CSDplot');xlim([1 2000]);
        end
        clear CSD LFP CSD_dum CSDplot LFPplot CSD_matr fl_write LFP_tbase sam_length
    end
    clear fl_length fl_name
end

   

    
            