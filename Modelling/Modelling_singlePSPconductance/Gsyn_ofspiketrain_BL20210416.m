%generates the postsynaptic gsyn timecourse as a response of a spike train
%sampling rate of the output is 20 kHz
clear all

Tsyn = 1.4;     % synaptic delay (ms)
Tdec = 35;      % tau decay (ms) 
Tris = 4.24;    % tau rise (ms)

%experiment data
bs_name='B'; bs_num='207'; bs_exp='b'; bs_typ='sp';                 % unit basic data entry
fl_ord=[1];                                                         % should be only one file

un_chns=[146];           %channel containing the spike times for the units
un_nums=[50];            %unit number (in kwik file)
un_shnk=[5];             %shanks from which unit was isolated
ch_LFP=[16];             %LFP cahannel in WAW.mat file

%reading data
fl_name=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord),'WAW.mat');         %file name
% reading spike times
ch_name=strcat(bs_name, bs_typ, bs_num, bs_exp,'_Ch', num2str(un_chns));
ch_load=load(num2str(fl_name),num2str(ch_name));
spks=(ch_load.(ch_name).('times')).*1000; clear ch_name ch_load;
% reading LFP channel
ch_name=strcat(bs_name, bs_typ, bs_num, bs_exp,'_Ch', num2str(ch_LFP));
ch_load=load(num2str(fl_name),num2str(ch_name));
fl_length =((((ch_load.(ch_name).('length'))-1)*(ch_load.(ch_name).('interval')))+(ch_load.(ch_name).('start')))*1000;
clear ch_name ch_load;

%generate a unitary postsynaptic conductance, detect the peak value for normalization
for je=1:10001
    tm=(je-1)/20;
    tmsc_u(je)=tm;
    gsyn_u(je)=(1-(exp((-1*tm)/Tris)))*(exp((-1*tm)/Tdec));
end; clear tm je
[pkval,pktim]=max(gsyn_u);
gsyn_u= gsyn_u./pkval;

%generate the gsyn base-trace from zeros at 20 kHz sampling rate
gsyn_0=zeros(1,(fl_length*20)+1);
gsyn_s=gsyn_0;
for je=1:length(spks)
    gsyn_a = gsyn_0;
    smpl_fst = round((spks(je) + Tsyn)*20);
    smpl_lst = smpl_fst + 10000;
    if smpl_lst<=size(gsyn_s,2)
        gsyn_a(1,smpl_fst:smpl_lst)=gsyn_u;
    end
    gsyn_s=gsyn_s+gsyn_a; 
    clear gsyn_a smpl_fst smpl_lst
end

%exporting datfile
%saving the relevant data
writename=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord),'_gsynPSC_s', num2str(un_shnk),'u',num2str(un_nums),'_Tsyn', num2str(Tsyn*100), '_Tdec', num2str(Tdec*100),'_Tris',num2str(Tris*100),'.dat');
fid=fopen(writename, 'w');
fwrite(fid,gsyn_s,'single');
fclose(fid);


