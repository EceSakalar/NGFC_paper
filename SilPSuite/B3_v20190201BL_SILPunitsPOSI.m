%This program reads original .dat files and takes waveform averages 
%for all channels aliogned based on the optimal channels (where the maximum amplitude is observed
%and estimates the position of the cell relative to the probe 
clear all; close all

%unit basic data entry
%textfile nameroot for shank 1 and shank 2
tx_name{1}='B193_d2_shL';
tx_name{2}='B193_d2_shC';
tx_name{3}='B193_d2_shR';

%datfile names for {shank#, file#}
dt_name{1,1}='B193b0_L.dat';
dt_name{2,1}='B193b0_C.dat';
dt_name{3,1}='B193b0_R.dat';

bs_name='B'; bs_num='193'; bs_exp='b'; bs_typ='sx'; bs_dayname='d02';

ch_n = [8, 16, 8];                                      %number of channels in the dat files
ch_rate = 20000;                                        %give the silicon probe sampling rate (in Hz)
%channels in the dat file with the maximum amplitude
ch_max = [6 6 5 7 6 7 5 7 4 1 4 6 3 4 4 3 6 3 3 14 11 11 11 11 10 4 7 4 6 7 3 3 6 2 4 6 7 3 6 4 3 6 6 7 7];
%channels in the spike2 files and originating tetrodes
ch_uns = [21:65];
%unit numbers (in kwik file)
un_nums=[19 20 30 40 51 58 74 111 113 114 115 116 118 122 123 124 125 126 128 62 68 75 84 94 95 20 27 32 70 85 91 98 140 142 147 148 149 158 161 163 164 165 168 169 170];                 
%shank identifier (should match the dat file order)
ch_tet =[repmat([1],1,19) repmat([2],1,6) repmat([3],1,20)];

%files to be included
fl_ord=[0];

%window for the spike extraction (in ms, one side), and resample rate for spike
wav_win = 3;
wav_rate = 100000;

tic
%unit stepper
for un_indx=1:length(ch_uns);
    
    %file stepper
    for fl_indx=1:length(fl_ord)
        
        %read the dat file
        acc_datname = dt_name{ch_tet(un_indx),fl_indx};
        dat_fid = fopen(acc_datname);
        dims = [ch_n(ch_tet(un_indx)),inf];
        acc_data = fread(dat_fid,dims,'int16=>int16'); acc_data=(double(acc_data')).*(10/65536);
        fclose(dat_fid);
        clear dat_fid dims acc_datname dims
       
        %read the spikes (for actual cell actual file)
        acc_txtname = strcat(tx_name{ch_tet(un_indx)}, '_u', num2str(un_nums(un_indx)), '_f', num2str(fl_ord(fl_indx)),'.txt');
        if dlmread(acc_txtname)==0;
            acc_spks=[];
        else
            acc_spks = dlmread(acc_txtname);
        end
        clear acc_txtname
        
        
        %collect & store all spikes for all units
        if fl_indx==1;
            spks{un_indx} = acc_spks;
            t_incr = size(acc_data,1)/ch_rate;
        else
            spks{un_indx} = [spks{un_indx}; (acc_spks+t_incr)];
            t_incr = t_incr+(size(acc_data,1)/ch_rate);
        end
        
        if un_indx == 1 && fl_indx == length(fl_ord)
            t_rec=t_incr;
        end
        
        
        %taking the actual channel for the spike extraction
        acc_chan = acc_data(:,(ch_max(un_indx)));
        %generating spike file for all the channels of the shank
        for chan=1:size(acc_data,2)
            sp_file{chan}=zeros(length(acc_spks),((wav_win*(ch_rate/1000))*2)+1);
        end
        
        if isempty(acc_spks)
            sp_file=[];
        else
            for sp_indx = 1:length(acc_spks)
                %extracting the wave on the original spike detected with 0.5 ms window
                acc_align = acc_chan((round(acc_spks(sp_indx)*ch_rate))-(0.5*(ch_rate/1000)):(round(acc_spks(sp_indx)*ch_rate))+(0.5*(ch_rate/1000)));
                %realign if not accurate
                [~,acc_minpos] = min(acc_align);
                acc_corr=acc_minpos-((0.5*(ch_rate/1000))+1); clear acc_align acc_minpos;
                %extract the new corrected spike
                if ((round(acc_spks(sp_indx)*ch_rate))-(wav_win*(ch_rate/1000))+acc_corr)>0&&((round(acc_spks(sp_indx)*ch_rate))+(wav_win*(ch_rate/1000))+acc_corr)<size(acc_chan,1)
                   acc_wave = acc_data((((round(acc_spks(sp_indx)*ch_rate))-(wav_win*(ch_rate/1000))+acc_corr):((round(acc_spks(sp_indx)*ch_rate))+(wav_win*(ch_rate/1000))+acc_corr)),:);
                    clear acc_corr
                    for chan=1:size(acc_data,2)
                        sp_file{chan}(sp_indx,:) = acc_wave(:,chan); 
                    end; clear acc_wave
                end
            end
        end
        clear sp_indx acc_spks
        
        %collect spike data over the files
        if fl_indx==1
            sp_all=sp_file;
        else
            for chan=1:size(acc_data,2)
            sp_all{chan}=[sp_all{chan}; sp_file{chan}];
            end
        end
        clear sp_file acc_chan acc_data
    end
    clear fl_indx t_incr
    
    char_units{un_indx,1} = ch_tet(un_indx);
    char_units{un_indx,2} = un_indx;
    char_units{un_indx,3} = sp_all;
    for chan=1:size(sp_all,2)
        char_units{un_indx,4}(:,chan) = mean(sp_all{chan});
    end
    
    %%ITT TARTOK
    %resampling & baseline correcting the spike (mean)
    char_units{un_indx,5}(1,:) = [(-1*wav_win):(1/wav_rate)*1000:wav_win];
    char_units{un_indx,5}(2,:) = interp1((-1*wav_win):(1/ch_rate)*1000:wav_win,mean(sp_all),(-1*wav_win):(1/wav_rate)*1000:wav_win)-mean(sp_all(:,1));
    char_units{un_indx,5}(3,:) = interp1((-1*wav_win):(1/ch_rate)*1000:wav_win,mean(sp_all),(-1*wav_win):(1/wav_rate)*1000:wav_win,'spline');
    %smoothing by ten samples, taking until 40%
    inter1 = (smooth(mean(sp_all),10));
    inter2 = inter1([[1:10:((round(length(mean(sp_all))*0.03333))*10)+1] floor(length(mean(sp_all))*0.8333)  floor(length(mean(sp_all))*0.916666) end]);
    inter3 = [(-1*wav_win):(1/ch_rate)*1000:wav_win];
    inter4 = inter3([[1:10:((round(length(mean(sp_all))*0.03333))*10)+1] floor(length(mean(sp_all))*0.8333)  floor(length(mean(sp_all))*0.916666) end]);
    char_units{un_indx,5}(4,:) = interp1(inter4,inter2,(-1*wav_win):(1/wav_rate)*1000:wav_win,'spline');
    char_units{un_indx,5}(5,:) = char_units{un_indx,5}(3,:) -  char_units{un_indx,5}(4,:); 
    clear sp_all inter1 inter2 inter3 inter4
    
    char_units{un_indx,6} = [t1 t2 aPeak sym];
    clear t1 t2 sym sPeak aPeak sFirst sLast 
    
end
clear un_indx
'first calculus finished'
toc



datum=date;

%save the relevant data
writename=strcat(bs_name, bs_num, bs_exp, bs_dayname, '_unitcharSILPR.mat');
save(writename, 'datum', 'ch_tet','ch_uns', 'un_nums','char_units','char_intact','char_intactQ','char_xcor','t_rec','fl_ord','spks');


