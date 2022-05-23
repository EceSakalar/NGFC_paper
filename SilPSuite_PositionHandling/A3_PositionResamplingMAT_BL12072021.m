clear all
%desired sampling interval
smpl_int=0.02;
%experiment data
bs_nam='B181g1';
ch_nam ='Bsx181g_Ch21';

%read the position time series
fl_name=strcat(bs_nam,'POS.mat');
ch_load=load(fl_name,ch_nam);

pos_orig=(ch_load.(ch_nam).values);
tim_orig=((ch_load.(ch_nam).start):(ch_load.(ch_nam).interval):(ch_load.(ch_nam).start)+((length(ch_load.(ch_nam).values)-1)*(ch_load.(ch_nam).interval)));
pts_orig=timeseries(pos_orig,tim_orig);
tim_new=(smpl_int:smpl_int:(smpl_int*(ceil(max(tim_orig)/smpl_int))));
pts_new=resample(pts_orig,tim_new);
pos_new=pts_new.Data;

exp_new=[[0 tim_new]' [pos_new(1,1); pos_new]];
if isnan(exp_new(end,2))
    exp_new(end,2)=exp_new(end-1,2);
end
    
%saving
fl_name=strcat(bs_nam, 'POSrs.txt');
dlmwrite(fl_name,exp_new,'precision',9);