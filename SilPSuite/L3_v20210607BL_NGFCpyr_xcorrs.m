%This program reads spike trains for given NGFCs and also pooled pyramidal cells
%and calculates cross corelograms to address significant interactions
%only periods of theta oscillations are included

clear all; %close all

%basic entry (only single file data is accepted)
bs.name='B'; bs.num='207'; bs.exp='b'; bs.typ='sp';
bs.flnum=[1];
bs.ngfname='s5u050';

%PYRAMIDAL CELLS
pyr.chns = [130,131,132,133,140,141,142,144,148,150,151,152,153];
%unit numbrs (in kwik file)
pyr.nums = [6,8,110,15,102,103,108,33,54,57,132,138,65];
%shank identifier (should match the dat file order)
pyr.shnk =[repmat([1],1,3) 2 repmat([3],1,3) repmat([4],1,1) repmat([5],1,4) repmat([6],1,1)];                       %shanks from which unit was isolated

%NGFCs
ngf.chn = [146];
%unit numbers

%Reading spike data
fl_name=strcat(bs.name, bs.num, bs.exp, num2str(bs.flnum),'WAW.mat'); %file name
%PYRAMIDAL CELLS
for pyr_indx=1:length(pyr.chns)
    ch_name=strcat(bs.name, bs.typ, bs.num, bs.exp,'_Ch', num2str(pyr.chns(pyr_indx)));
    ch_load=load(num2str(fl_name),num2str(ch_name));
    spks_pyr{pyr_indx}=(ch_load.(ch_name).('times')); clear ch_name ch_load;
    if pyr_indx==1
        spks_pyrpool=spks_pyr{pyr_indx};
    else
        spks_pyrpool=[spks_pyrpool;spks_pyr{pyr_indx}];
    end
end; clear pyr_indx

%NGFCs
ch_name=strcat(bs.name, bs.typ, bs.num, bs.exp,'_Ch', num2str(ngf.chn));
ch_load=load(num2str(fl_name),num2str(ch_name));
spks_ngf=(ch_load.(ch_name).('times')); clear ch_name ch_load;
clear fl_name;

%reading theta periods
th_pname=strcat(bs.name, bs.num, bs.exp, num2str(bs.flnum),'THP.txt');%theta periods
th_pers=dlmread(th_pname); clear th_pname;

%Limiting spikes for the time span of NGFC recording
spks_pyrngf=spks_pyrpool(spks_pyrpool>=(min(spks_ngf)-0.5)&spks_pyrpool<=(max(spks_ngf)+0.5));

%limiting spikes to theta periods
spks_pyrtht=[];spks_ngftht=[];
for ka=1:size(th_pers,1)  %theta period stepper
    spks_Nact = spks_ngf(spks_ngf>=th_pers(ka,1) & spks_ngf<=th_pers(ka,2));
    spks_Pact = spks_pyrngf(spks_pyrngf>=th_pers(ka,1) & spks_pyrngf<=th_pers(ka,2));
    if isempty(spks_ngftht)
        spks_ngftht=spks_Nact;
    elseif ~isempty(spks_Nact)
        spks_ngftht=[spks_ngftht; spks_Nact];
    end
    if isempty(spks_pyrtht)
        spks_pyrtht=spks_Pact;
    elseif ~isempty(spks_Pact)
        spks_pyrtht=[spks_pyrtht; spks_Pact];
    end
    clear spks_Nact spks_Pact
end; clear ka
spks_pyrtht=sort(spks_pyrtht);

%do the NGFC(presynaptic) PYR(postsynaptic) cross-correlograms
%if presynaptic partner (NGF) is n<10000 do it straight
if length(spks_ngftht)<10000
    preSPKnum = length(spks_ngftht);
    pstSPKnum = length(spks_pyrtht);
    x_mat = repmat(spks_pyrtht,1,size(spks_ngftht,1)) - repmat(spks_ngftht',size(spks_pyrtht,1),1);
    x_mat(abs(x_mat)>1)=NaN;
    [~,V_indx]=find(~isnan(x_mat));
    M_indx=find(~isnan(x_mat));
    x_vect=x_mat(M_indx);
    clear M_indx x_mat
    %if presyn partner is big (>10000) segment it!
else
    seg_num=ceil(length(spks_ngftht)/9000);
    preSPKnum = length(spks_ngftht);
    pstSPKnum = length(spks_pyrtht);
    x_vect = NaN;V_indx = NaN;
    for seg_indx=1:seg_num
        if seg_indx==1
            seg_spks = spks_ngftht((((seg_indx-1)*9000)+1):(seg_indx*9000));
            seg_mat = repmat(spks_pyrtht,1,size(seg_spks,1)) - repmat(seg_spks',size(spks_pyrtht,1),1);
            seg_mat(abs(seg_mat)>1)=NaN;
            [~,V_temp] = find(~isnan(seg_mat));
            if ~isempty(V_temp)
                V_indx = [V_indx; V_temp];
            end
            clear V_temp
            M_indx = find(~isnan(seg_mat));
            x_temp = seg_mat(M_indx);
            clear M_indx
            if ~isempty(x_temp)
                x_vect = [x_vect; x_temp];
            end
            clear x_temp seg_mat seg_spks
            
        elseif seg_indx==seg_num
            seg_spks=spks_ngftht((((seg_indx-1)*9000)+1):end);
            seg_mat = repmat(spks_pyrtht,1,size(seg_spks,1)) - repmat(seg_spks',size(spks_pyrtht,1),1);
            seg_mat(abs(seg_mat)>1)=NaN;
            [~,V_temp] = find(~isnan(seg_mat));
            V_temp = V_temp + (9000*(seg_indx-1));
            if ~isempty(V_temp)
                V_indx = [V_indx; V_temp];
            end
            clear V_temp
            M_indx = find(~isnan(seg_mat));
            x_temp = seg_mat(M_indx);
            clear M_indx
            if ~isempty(x_temp)
                x_vect = [x_vect; x_temp];
            end
            clear x_temp seg_mat seg_spks
            
        else
            seg_spks=spks_ngftht((((seg_indx-1)*9000)+1):(seg_indx*9000));
            seg_mat = repmat(spks_pyrtht,1,size(seg_spks,1)) - repmat(seg_spks',size(spks_pyrtht,1),1);
            seg_mat(abs(seg_mat)>1)=NaN;
            [~,V_temp] = find(~isnan(seg_mat));
            V_temp = V_temp + (4000*(seg_indx-1));
            if ~isempty(V_temp)
                V_indx = [V_indx; V_temp];
            end
            clear V_temp
            M_indx = find(~isnan(seg_mat));
            x_temp = seg_mat(M_indx);
            clear M_indx
            if ~isempty(x_temp)
                x_vect = [x_vect; x_temp];
            end
            clear x_temp seg_mat
        end
    end; clear seg_indx seg_num seg_spks
    x_vect(1)=[];
    V_indx(1)=[];
end

%filling the xcor
x_cor = (((hist(x_vect,-0.101:0.001:0.101))./((preSPKnum))).*1000)';
x_cor(1)=[]; x_cor(end)=[];

%generate shuffling shu_xcor variable
shu_xcor=NaN(100,203);
for je=1:100
    jit_nums=((rand(1,preSPKnum)).*0.02)-0.01;
    jit_vect=x_vect+jit_nums(V_indx)';
    shu_xcor(je,:) = hist(jit_vect,-0.101:0.001:0.101);
    clear jit_nums jit_vect
end; clear je
shu_xcor(:,1)=[]; shu_xcor(:,end)=[];
clear V_indx x_vect

%calculus in the forward direction
%forward direction--------------------------------------

x_cor(:,2) = ((mean(shu_xcor,1))./(size(spks_ngftht,1))).*1000;
x_cor(:,3) = ((std(shu_xcor,1,1))./(size(spks_ngftht,1))).*1000;
x_cor(:,4) = x_cor(:,2)+(3.3.*(x_cor(:,3)));
x_cor(:,5) = x_cor(:,2)-(3.3.*(x_cor(:,3)));
x_cor(:,6) = NaN(201,1);                                        %over or under confidence
x_cor(x_cor(:,1)>x_cor(:,4),6)=1;
x_cor(x_cor(:,1)<x_cor(:,5),6)=-1;
x_cor(:,7) = (((x_cor(:,1))./(x_cor(:,2)))-1)*100;              %excess or missing (percent of expected)
x_cor(:,8) = (((x_cor(:,1))-(x_cor(:,2))).*preSPKnum)./1000;    %excess or missing (spikes)

%%--------------------------PLOTTING the data
figure
%textual

subplot (3,1,1)
axis off
text(0,9,strcat(bs.name, bs.num, bs.exp,'ONLY THETA'));
text(0,6,strcat('Pre unit: NGFC', bs.ngfname));
text(0,3,strcat('Post units pooled PYRs, n=', num2str(length(pyr.chns))));
text(0,0,strcat('100 x shuffle, jitter +/- 10 ms'))
axis ([0 1 0 12])

%crosscorelations, expected rates and confidence intervals from shuffling
subplot (3,1,2)
bar(-100:100,x_cor(1:201,1))
hold on
line(-100:100,x_cor(1:201,2),'Color',[1 0 0]);
line(-100:100,x_cor(1:201,4),'Color',[0.5 0.5 0.5]);
line(-100:100,x_cor(1:201,5),'Color',[0.5 0.5 0.5]);
line([0 0],[0 1.1*(max(max(x_cor(1:201,1:5))))],'Color',[0 1 0]);
line(-100:100,abs(x_cor(1:201,6))*(1.05*(max(max(x_cor(1:201,1:5))))),'Color',[0 0 0],'Marker','.') 
xlabel ('time from NGFC spike (ms)'); ylabel('Pyr pooled firing (Hz)');
axis([-100, 100, 0, 1.1*(max(max(x_cor(1:201,1:5))))]);
hold off

%percentages (red) and spikes (blue)
subplot (3,1,3)
hold on
line(-100:100,x_cor(1:201,7),'Color',[1 0 0]);
line(-100:100,x_cor(1:201,8),'Color',[0 0 1]);
line([-100 100],[30 30],'Color',[0.5 0.5 0.5]);
line([-100 100],[-30 -30],'Color',[0.5 0.5 0.5]);
xlabel ('time after NGFC spike (ms)'); ylabel('excess/missing %/# spikes');
axis([-100, 100, (-20+(min(min(x_cor(1:201,7:8))))), (20+(max(max(x_cor(91:201,7:8)))))]);
hold off

datum=date;

%save the relevant data
writename=strcat(bs.name, bs.num, bs.exp, '_NGFC_', bs.ngfname,'vs',num2str(length(pyr.chns)),'PYRs.mat');
save(writename, 'datum', 'bs','ngf', 'pyr','x_cor','shu_xcor','-v7.3');

