clear all
%generate firing rate measurements for Theta, SWR, & non-theta-non-swr periods
fl_name='B187b1WAW.mat';
ch_name='Bsx187b_Ch43';
bs_name='B187b1';

%loading the spike train
spk_load=load(fl_name,ch_name);
spks=(spk_load.(ch_name).times);

%THETA
t_name=strcat(bs_name,'THP.txt');
t_pers=dlmread(t_name);
%time
t_time=(t_pers(:,2)-t_pers(:,1));
%spike number
for je=1:size(t_pers,1)
    t_spks(je,1)=length(spks(spks>=t_pers(je,1)&spks<t_pers(je,2)));
end
t_rate=t_spks./t_time;

%SWR
s_name=strcat(bs_name,'SWP.txt');
s_pers=dlmread(s_name);
%time
s_time=(s_pers(:,2)-s_pers(:,1));
%spike number
for je=1:size(s_pers,1)
    s_spks(je,1)=length(spks(spks>=s_pers(je,1)&spks<s_pers(je,2)));
end
s_rate=s_spks./s_time;

%non-THETA-non-SWR
n_name=strcat(bs_name,'NTSP.txt');
n_pers=dlmread(n_name);
%time
n_time=(n_pers(:,2)-n_pers(:,1));
%spike number
for je=1:size(n_pers,1)
    n_spks(je,1)=length(spks(spks>=n_pers(je,1)&spks<n_pers(je,2)));
end
n_rate=n_spks./n_time;

stats(1,1)=sum(t_time);
stats(2,1)=sum(s_time);
stats(3,1)=sum(n_time);

stats(5,1)=sum(t_spks);
stats(6,1)=sum(s_spks);
stats(7,1)=sum(n_spks);

stats(9,1)=stats(5,1)/stats(1,1);
stats(10,1)=stats(6,1)/stats(2,1);
stats(11,1)=stats(7,1)/stats(3,1);


