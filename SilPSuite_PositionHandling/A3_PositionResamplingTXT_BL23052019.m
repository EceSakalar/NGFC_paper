%desired sampling interval
smpl_int=0.02;
%experiment data
bs_nam='B'; bs_num='204'; bs_exp='a1';

%resampling
fl_name=strcat(bs_nam, bs_num, bs_exp, 'POS.txt');
pos_orig=dlmread (fl_name);
pos_orig_ts=timeseries(pos_orig(:,3),pos_orig(:,1));
pos_rsmp_tm=[(smpl_int*(ceil(min(pos_orig_ts.Time)/smpl_int))):smpl_int:(smpl_int*(floor(max(pos_orig_ts.Time)/smpl_int)))];
pos_rsmp_ts=resample(pos_orig_ts,pos_rsmp_tm);
pos_rsmp=[pos_rsmp_ts.Time pos_rsmp_ts.Data]; 

%saving
fl_name=strcat(bs_nam, bs_num, bs_exp, 'POSrs.txt');
dlmwrite(fl_name,pos_rsmp);