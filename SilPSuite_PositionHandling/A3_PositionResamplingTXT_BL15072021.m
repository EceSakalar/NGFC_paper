clear all
%desired sampling interval
smpl_int=0.02;
%experiment data
bs_nam='B204b1';

%read the position time series
fl_name=strcat(bs_nam, 'POS.txt');
pos_read=dlmread(fl_name);

pos_orig=pos_read(:,2);
tim_orig=pos_read(:,1);
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