%This script works if all data of the experiment is in a single file!!!
%CALCULATES
%-characteristics of the spikes recorded by the glass electrodes (only one cell)
%-Cross correlograms with other cells recorded from the silicon probe for theta periods only!!

clear all; %close all

%unit basic data entry
bs_name='B'; bs_num='186'; bs_exp='d'; bs_typ='sx'; bs_dayname='d4';
%the order of files
fl_ord=[1];
%Give the channels for the unit spikes for the glass cell & the waveform
gl_evnt = 55;   
gl_wave = 33;

%channels in the spike2 files
ch_uns = [35:53];
%unit numbers (in kwik file)
un_nums=[0,5,7,19,23,77,85,104,170,172,202,209,211,50,56,58,73,195,196];                 
%shank identifier (should match the dat file order)
un_shnk =[repmat([1],1,13) repmat([2],1,2) repmat([3],1,4)];                       %shanks from which unit was isolated

%parameters for extracting the waveforms
%window for the spike extraction (in ms, one side), and resample rate for spike
wav_win = 4;
wav_rate = 100000;


tic
fl_name=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord),'WAW.mat');   %file name
th_pname=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord),'THP.txt');  %theta periods
th_pers=dlmread(th_pname); clear th_pname;

%---------READING SPIKES
%read the spikes for the glass unit
ch_name = strcat(bs_name, bs_typ, bs_num, bs_exp,'_Ch', num2str(gl_evnt));
ch_load=load(num2str(fl_name),num2str(ch_name));
glassp_all = (ch_load.(ch_name).('times')); clear ch_name ch_load;
glassp_tht=[];
for ka=1:size(th_pers,1)  %theta period stepper
    spks_act = glassp_all(glassp_all>=th_pers(ka,1) & glassp_all<=th_pers(ka,2));
    if isempty(glassp_tht)
        glassp_tht=spks_act;
    elseif ~isempty(spks_act)
        glassp_tht=[glassp_tht; spks_act];
    end
    clear spks_act
end; clear ka

%read the spikes for interaction partners
for un_indx=1:length(ch_uns)
    ch_name = strcat(bs_name, bs_typ, bs_num, bs_exp,'_Ch', num2str(ch_uns(un_indx)));
    ch_load=load(num2str(fl_name),num2str(ch_name));
    %read the spikes (for actual cell actual file)
    spks{un_indx} = (ch_load.(ch_name).('times')); clear ch_name ch_load;
    
    spks_th{un_indx}=[];
    for ka=1:size(th_pers,1)  %theta period stepper
        spks_act = spks{un_indx}(spks{un_indx}>=th_pers(ka,1) & spks{un_indx}<=th_pers(ka,2));
        if isempty(spks_th{un_indx})
            spks_th{un_indx}=spks_act;
        elseif ~isempty(spks_act)
            spks_th{un_indx}=[spks_th{un_indx}; spks_act];
        end
        clear spks_act
    end; clear ka
end; clear un_indx

%%%%%----------GLASS UNIT characteristics

%read the mat file for spike voltage recording
ch_name=strcat(bs_name, bs_typ, bs_num, bs_exp,'_Ch', num2str(gl_wave));
ch_load=load(num2str(fl_name),num2str(ch_name));
gl_chan=(ch_load.(ch_name).('values'));
gl_ofst=(ch_load.(ch_name).('offset'));
gl_ival=(ch_load.(ch_name).('interval')); clear ch_load ch_name;
t_rec=max(glassp_all)-min(glassp_all);

%collect and store spike shapes
gl_wfrm=gl_chan(round((repmat(round(glassp_all./gl_ival),1,round(((wav_win*(0.001/gl_ival))*2)+1)))+(repmat([((-1*wav_win)/(gl_ival*1000)):1:((1*wav_win)/(gl_ival*1000))],length(glassp_all),1))));
char_GLunits{1} = 'glass';
char_GLunits{2} = gl_wfrm;
char_GLunits{3} = mean(gl_wfrm);
    
%resampling & baseline correcting the spike (mean)
char_GLunits{5}(1,:) = [(-1*wav_win):(1/wav_rate)*1000:wav_win];
char_GLunits{5}(2,:) = interp1((-1*wav_win):(gl_ival*1000):wav_win,mean(gl_wfrm),(-1*wav_win):(1/wav_rate)*1000:wav_win)-mean(gl_wfrm(:,1));
char_GLunits{5}(3,:) = interp1((-1*wav_win):(gl_ival*1000):wav_win,mean(gl_wfrm),(-1*wav_win):(1/wav_rate)*1000:wav_win,'spline');

%smoothing by ten samples, taking until 40%
inter1 = (smooth(mean(gl_wfrm),10));
inter2 = inter1([[1:10:((round(length(mean(gl_wfrm))*0.03333))*10)+1] floor(length(mean(gl_wfrm))*0.8333)  floor(length(mean(gl_wfrm))*0.916666) end]);
inter3 = [(-1*wav_win):gl_ival*1000:wav_win];
inter4 = inter3([[1:10:((round(length(mean(gl_wfrm))*0.03333))*10)+1] floor(length(mean(gl_wfrm))*0.8333)  floor(length(mean(gl_wfrm))*0.916666) end]);
char_GLunits{5}(4,:) = interp1(inter4,inter2,(-1*wav_win):(1/wav_rate)*1000:wav_win,'spline');
char_GLunits{5}(5,:) = char_GLunits{5}(3,:) -  char_GLunits{5}(4,:);
clear sp_all inter1 inter2 inter3 inter4

%spike charistics calculus
%find the peak and the trough
[aPeak,sPeak] = max(char_GLunits{5}(5,:));
[aTrgh,sTrgh] = min(char_GLunits{5}(5,:));

%find 'decay' to 0.9 before the peak and after the trough of spike
sFirst = find(char_GLunits{5}(5,1:sPeak)<=aPeak*0.1, 1, 'last' );
sLast = find(char_GLunits{5}(5,sTrgh:end)>=aTrgh*0.1, 1, 'first' );
t1 = ((sFirst-sPeak)/wav_rate)*1000;
t2 = (sLast/wav_rate)*1000;
sym =(abs(t1)-abs(t2))/(abs(t1)+abs(t2));
t_all = t2-t1;

char_GLunits{6} = [t1 t2 t_all aPeak sym];
clear t1 t2 t_all sym sPeak aPeak sFirst sLast sTrgh aTrgh
clear gl_indx t_incr
'first calculus finished'
toc
   
%generate autocorrelograms
%make autocorrelogram via automatrix (+/- 100 ms, 1ms bins)
%make an accomodating matrix
a_mat=NaN(length(glassp_all),201);
%do it straight
if length(glassp_all)<5000
    seg_mat0 = [zeros(100,size(glassp_all,1)); repmat(glassp_all,1,size(glassp_all,1))-repmat(glassp_all',size(glassp_all,1),1); zeros(100,size(glassp_all,1))];
    for ka=0:-1:-200
        a_mat(:,((ka*-1)+1))=diag(seg_mat0,ka);
    end; clear ka seg_mat0
    
    %do it chunked
else
    seg_num=ceil(length(glassp_all)/4000);
    for seg_indx=1:seg_num
        if seg_indx==1
            seg_spks=glassp_all((((seg_indx-1)*4000)+1):(seg_indx*4000));
            seg_mat1=repmat(seg_spks,1,size(seg_spks,1))-repmat(seg_spks',size(seg_spks,1),1);
            seg_mat2=[zeros(100,4000); seg_mat1; zeros(100,4000)];
            seg_mat3=NaN(4000,201);
            for ka=0:-1:-200
                seg_mat3(:,((ka*-1)+1))=diag(seg_mat2,ka);
            end; clear ka seg_spks seg_mat1 seg_mat2
            a_mat(1:3900,:)=seg_mat3(1:3900,:);
        elseif seg_indx==seg_num
            seg_spks=glassp_all((((seg_indx-1)*4000)-99):end);
            seg_mat1=repmat(seg_spks,1,size(seg_spks,1))-repmat(seg_spks',size(seg_spks,1),1);
            seg_mat2=[zeros(100,size(seg_mat1,2)); seg_mat1; zeros(100,size(seg_mat1,2))];
            seg_mat4=NaN(size(seg_mat1,1),201);
            for ka=0:-1:-200
                seg_mat4(:,((ka*-1)+1))=diag(seg_mat2,ka);
            end; clear ka seg_spks seg_mat1 seg_mat2
            seg_matb=seg_mat3((size(seg_mat3,1)-99):end,:);
            [indic1, indic2]=find(seg_matb==0);
            seg_matb(indic1,indic2)=seg_mat4(indic1,indic2);
            seg_mat4(1:100,:)=seg_matb; clear seg_matb indic1 indic2
            a_mat((((seg_indx-1)*4000)-99):end,:)=seg_mat4(1:end,:);
            clear seg_mat3 seg_mat4
        else
            seg_spks=glassp_all((((seg_indx-1)*4000)-99):(seg_indx*4000));
            seg_mat1=repmat(seg_spks,1,size(seg_spks,1))-repmat(seg_spks',size(seg_spks,1),1);
            seg_mat2=[zeros(100,4100); seg_mat1; zeros(100,4100)];
            seg_mat4=NaN(4100,201);
            for ka=0:-1:-200
                seg_mat4(:,((ka*-1)+1))=diag(seg_mat2,ka);
            end; clear ka seg_spks seg_mat1 seg_mat2
            seg_matb=seg_mat3((size(seg_mat3,1)-99):end,:);
            [indic1, indic2]=find(seg_matb==0);
            seg_matb(indic1,indic2)=seg_mat4(indic1,indic2);
            seg_mat4(1:100,:)=seg_matb; clear seg_matb indic1 indic2
            a_mat(((((seg_indx-1)*4000)-99):(seg_indx*4000)-100),:)=seg_mat4(1:4000,:);
            seg_mat3=seg_mat4; clear seg_mat4
        end
    end; clear seg_indx seg_num
end
a_cor = (sum(hist(a_mat,-0.101:0.001:0.101),2))./(size(glassp_all,1)).*1000;
a_cor(1)=[]; a_cor(end)=[];a_cor(101)=0;
%for the second column take only values >half maximum
a_corS = a_cor; a_corS(a_cor<(max(a_cor)/2))=NaN;
char_GLunits{7} = [a_cor a_corS];
%calculate the mean rate, n and weighted avergae of >mean(rate) time values
char_GLunits{8} = [length(glassp_all)/t_rec length(glassp_all) (nansum(abs([-100:1:100]').*a_corS))/nansum(a_corS)];
clear a_mat a_cor a_corS
'autocorrelograms done'
toc

%------Plot unit spike and autocorrelogram results
figure
%unit identity
subplot (3,1,1)
axis off
text(0,7,strcat(bs_name, bs_num, bs_exp,'(',bs_dayname,')'));
text(0,5,strcat('unit:glass in ch', num2str(gl_evnt)));
axis ([0 1 0 8])
%spike and spike parameters
subplot (3,1,2)
plot (char_GLunits{5}(1,:),char_GLunits{5}(2:5,:))
set (gca,'ColorOrder',[0.2 0.2 0.2; 0.2 0.2 0.2;0.2 0.2 0.2;1 1 1]);
axis tight
xlabel ('time (ms)');ylabel('amplitude');
text (-2.9 , (max(max(char_GLunits{5}(2:5,:)))*0.9) , strcat('t1=',num2str(char_GLunits{6}(1))));
text (-2.9 , (max(max(char_GLunits{5}(2:5,:)))*0.6) , strcat('t2=',num2str(char_GLunits{6}(2))));
text (-2.9 , (max(max(char_GLunits{5}(2:5,:)))*0.3) , strcat('sym=',num2str(char_GLunits{6}(5))));
text (-2.9 , (max(max(char_GLunits{5}(2:5,:)))*0.0) , strcat('t_all=',num2str(char_GLunits{6}(3))));
%expand the graph
pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi)
clear pozi

%autocorrelogram and parameters
subplot (3,1,3)
bar(-100:1:100,char_GLunits{7});
xlabel ('time (ms)');ylabel('rate (Hz)');
axis([-50, 50, 0, max(char_GLunits{7}(:,1))]);
text (-49 , (max(char_GLunits{7}(:,1))*0.95) , strcat('rate=',num2str(char_GLunits{8}(1))));
text (-49 , (max(char_GLunits{7}(:,1))*0.75) , strcat('n=',num2str(char_GLunits{8}(2))));
text (-49 , (max(char_GLunits{7}(:,1))*0.55) , strcat('hub=',num2str(char_GLunits{8}(3))));
%expand the graph
pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi)
clear pozi

%---------CROSS CORRELOGRAMS
%make cross correlograms via cross matrices
%loop for the silicon probe partner
for pt_indx = 1:length(ch_uns)
    %if presynaptic partner is n<5000 do it straight
    if length(glassp_tht)<5000
        preSPKnum = length(glassp_tht);
        pstSPKnum = length(spks_th{pt_indx});
        x_mat = repmat(spks_th{pt_indx},1,size(glassp_tht,1)) - repmat(glassp_tht',size(spks_th{pt_indx},1),1);
        x_mat(abs(x_mat)>1)=NaN;
        [~,V_indx]=find(~isnan(x_mat));
        M_indx=find(~isnan(x_mat));
        x_vect=x_mat(M_indx);
        clear M_indx x_mat
        
        %if presyn partner is big segment it!
    else
        seg_num=ceil(length(glassp_tht)/4000);
        preSPKnum = length(glassp_tht);
        pstSPKnum = length(spks_th{pt_indx});
        x_vect = NaN;V_indx = NaN;
        for seg_indx=1:seg_num
            if seg_indx==1
                seg_spks = glassp_tht((((seg_indx-1)*4000)+1):(seg_indx*4000));
                seg_mat = repmat(spks_th{pt_indx},1,size(seg_spks,1)) - repmat(seg_spks',size(spks_th{pt_indx},1),1);
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
                seg_spks=glassp_tht((((seg_indx-1)*4000)+1):end);
                seg_mat = repmat(spks_th{pt_indx},1,size(seg_spks,1)) - repmat(seg_spks',size(spks_th{pt_indx},1),1);
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
                clear x_temp seg_mat seg_spks
                
            else
                seg_spks=glassp_tht((((seg_indx-1)*4000)+1):(seg_indx*4000));
                seg_mat = repmat(spks_th{pt_indx},1,size(seg_spks,1)) - repmat(seg_spks',size(spks_th{pt_indx},1),1);
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
    x_cor = ((hist(x_vect,-0.101:0.001:0.101))./(preSPKnum).*1000)';
    x_cor(1)=[]; x_cor(end)=[];
    char_xcorFORW{pt_indx}=x_cor;
    
    %filling the 'mirror' x_cor
    x_cor = ((hist((x_vect.*-1),-0.101:0.001:0.101))./(pstSPKnum).*1000)';
    x_cor(1)=[]; x_cor(end)=[];
    char_xcorBACK{pt_indx}=x_cor;
    clear x_cor
    
    %shuffling calculus for xcor---------------------------------------------
    %FIRST layer of INTERACTION matrix
    %NaN the mean number of spikes in all bins is <= 1
    %0 if there is no interaction
    %1 if there is significant interaction
    %SECOND layer of INTERACTION matrix
    %1 if COFIRING
    %THIRD layer of INTERACTION matrix
    %1 if EXCITATION
    %FOURTH layer of INTERACTION matrix
    %1 if INHIBITION
    
    shu_xcor=NaN(100,203);
    for je=1:100
        jit_nums=((rand(1,preSPKnum)).*0.02)-0.01;
        jit_vect=x_vect+jit_nums(V_indx)';
        shu_xcor(je,:) = hist(jit_vect,-0.101:0.001:0.101);
        clear jit_nums jit_vect
    end; clear je
    shu_xcor(:,1)=[]; shu_xcor(:,end)=[];
    toc
    clear V_indx x_vect
    
    %forward direction--------------------------------------
    x_cor = char_xcorFORW{pt_indx};
    x_cor(:,2) = (mean(shu_xcor,1)./(size(glassp_tht,1))).*1000;
    x_cor(:,3) = (std(shu_xcor,1,1)./(size(glassp_tht,1))).*1000;
    x_cor(:,4) = x_cor(:,2)+(3.3.*(x_cor(:,3)));
    x_cor(:,5) = x_cor(:,2)-(3.3.*(x_cor(:,3)));
    x_cor(:,6) = NaN(201,1);                                        %over or under confidence
    x_cor(x_cor(:,1)>x_cor(:,4),6)=1;
    x_cor(x_cor(:,1)<x_cor(:,5),6)=-1;
    x_cor(:,7) = (((x_cor(:,1))./(x_cor(:,2)))-1)*100;              %excess or missing (percent of expected)
    x_cor(:,8) = (((x_cor(:,1))-(x_cor(:,2))).*preSPKnum)./1000;    %excess or missing (spikes)
    char_xcorFORW{pt_indx}=x_cor;
    %char_xcorjit{un_indx,pt_indx}=shu_xcor;
    
    %check if any significant bins in the forward direction(for pre - post interaction)
    %for significant interactions store the quantification in INTACTQ matrix layers 1 (cofir), 2 (excitation, and 3 (inhibition)
    if nansum(abs(x_cor(101:end,6)))>0
        char_intact{1,pt_indx,1} = 1;
        char_intact{1,pt_indx,2} = 0;
        char_intact{1,pt_indx,3} = 0;
        char_intact{1,pt_indx,4} = 0;
        %check cofiring
        if x_cor(101,1)>x_cor(101,4) && x_cor(101,7)>30 && x_cor(101,8)>30;
            char_intact{1,pt_indx,2} = 1;
            char_intactQ{1,pt_indx,1} = (x_cor(101,8))/preSPKnum;
        end
        x_eval=x_cor(102:105,:);
        %check excitation
        if ~isempty(x_eval(x_eval(:,6)==1,6))
            if max(x_eval(x_eval(:,6)==1,7))>=30 && sum(x_eval(x_eval(:,6)==1,8))>=30;
                char_intact{1,pt_indx,3} = 1;
                char_intactQ{1,pt_indx,2} = sum(x_eval(x_eval(:,6)==1,8))/preSPKnum;
            end
        end
        %check inhibition
        if ~isempty(x_eval(x_eval(:,6)==-1,6))
            if min(x_eval(x_eval(:,6)==-1,7))<=-30 && sum(x_eval(x_eval(:,6)==-1,8))<=-30;
                char_intact{1,pt_indx,4} = 1;
                char_intactQ{1,pt_indx,3} = sum(x_eval(x_eval(:,6)==-1,8))/preSPKnum;
            end
        end
        clear x_eval
        strcat('pairs pre ', num2str(un_nums(pt_indx)), 'to post glass unit ')
        toc
    else
        char_intact{1,pt_indx,1} = 0;
        char_intact{1,pt_indx,2} = 0;
        char_intact{1,pt_indx,3} = 0;
        char_intact{1,pt_indx,4} = 0;
    end
    
    %inverted direction--------------------------------------
    x_cor = char_xcorBACK{pt_indx};
    temp1 = (mean(shu_xcor,1)./(size(spks_th{pt_indx},1))).*1000;
    temp2 = temp1;
    for ka=1:length(temp1)
        temp2(ka) = temp1(length(temp1)-ka+1);
    end
    x_cor(:,2) = temp2; clear temp1 temp2 ka
    
    temp1 = (std(shu_xcor,1,1)./(size(spks_th{pt_indx},1))).*1000;
    temp2 = temp1;
    for ka=1:length(temp1)
        temp2(ka) = temp1(length(temp1)-ka+1);
    end
    x_cor(:,3) = temp2; clear temp1 temp2 ka
    x_cor(:,4) = x_cor(:,2)+(3.3.*(x_cor(:,3)));
    x_cor(:,5) = x_cor(:,2)-(3.3.*(x_cor(:,3)));
    x_cor(:,6) = NaN(201,1);                                        %over or under confidence
    x_cor(x_cor(:,1)>x_cor(:,4),6)=1;
    x_cor(x_cor(:,1)<x_cor(:,5),6)=-1;
    x_cor(:,7) = (((x_cor(:,1))./(x_cor(:,2)))-1)*100;              %excess or missing (percent of expected)
    x_cor(:,8) = (((x_cor(:,1))-(x_cor(:,2))).*pstSPKnum)./1000;    %excess or missing (spikes)
    char_xcorBACK{pt_indx}=x_cor;
    
    %check if any significant bins in the reverse direction(for pre - post interaction)
    %for significant interactions store the quantification in INTACTQ matrix layers 1 (cofir), 2 (excitation, and 3 (inhibition)
    if nansum(abs(x_cor(101:end,6)))>0
        char_intact{2,pt_indx,1} = 1;
        char_intact{2,pt_indx,2} = 0;
        char_intact{2,pt_indx,3} = 0;
        char_intact{2,pt_indx,4} = 0;
        %check cofiring
        if x_cor(101,1)>x_cor(101,4) && x_cor(101,7)>30 && x_cor(101,8)>30;
            char_intact{2,pt_indx,2} = 1;
            char_intactQ{2,pt_indx,1} = (x_cor(101,8))/pstSPKnum;
        end
        x_eval=x_cor(102:105,:);
        %check excitation
        if ~isempty(x_eval(x_eval(:,6)==1,6))
            if max(x_eval(x_eval(:,6)==1,7))>=30 && sum(x_eval(x_eval(:,6)==1,8))>=30;
                char_intact{2,pt_indx,3} = 1;
                char_intactQ{2,pt_indx,2} = sum(x_eval(x_eval(:,6)==1,8))/pstSPKnum;
            end
        end
        %check inhibition
        if ~isempty(x_eval(x_eval(:,6)==-1,6))
            if min(x_eval(x_eval(:,6)==-1,7))<=-30 && sum(x_eval(x_eval(:,6)==-1,8))<=-30;
                char_intact{2,pt_indx,4} = 1;
                char_intactQ{2,pt_indx,3} = sum(x_eval(x_eval(:,6)==-1,8))/pstSPKnum;
            end
        end
        clear x_eval
        strcat('pairs pre glass unit to post ', num2str(un_nums(pt_indx)))
        toc
    else
        char_intact{2,pt_indx,1} = 0;
        char_intact{2,pt_indx,2} = 0;
        char_intact{2,pt_indx,3} = 0;
        char_intact{2,pt_indx,4} = 0;
    end
    clear x_cor shu_xcor
    
end; clear pt_indx
clear  preSPKnum pstSPKnum 


%------------PLOT CROSSCORELOGRAMS
%plot index
fg_indx=1; pl_indx=1;
for pr_indx=1:2       %glass unit loop pre or post
    %loop for the xcorr partner
    for pt_indx = 1:length(ch_uns)
        
        %new figure?
        if pl_indx==1
            figure
        end
        
        %units involved identity, interactions
        subplot (5,3,(((pl_indx-1)*3)+1))
        axis off
        text(0,11,strcat(bs_name, bs_num, bs_exp,'(',bs_dayname,')', ' ONLY THETA'));
        if pr_indx==1
            text(0,9,'Pre unit: glass');
            text(0,7,strcat('Post unit: shank',num2str(un_shnk(pt_indx)),' unit', num2str(un_nums(pt_indx))));
        elseif pr_indx==2
            text(0,9,strcat('Pre unit: shank',num2str(un_shnk(pt_indx)),' unit', num2str(un_nums(pt_indx))));
            text(0,7,'Post unit: glass');
        end
    
        text(0,5,strcat('Cofiring:',num2str(char_intact{pr_indx,pt_indx,2})));
        text(0,3,strcat('Excitation:',num2str(char_intact{pr_indx,pt_indx,3})));
        text(0,1,strcat('Inhibition:',num2str(char_intact{pr_indx,pt_indx,4})));
        axis ([0 1 0 12])
        
        %crosscorelations, expected rates and confidence intervals from shuffling
        if pr_indx==1
            x_cor=char_xcorFORW{pt_indx};
        elseif pr_indx==2
            x_cor=char_xcorBACK{pt_indx};
        end
                
        subplot (5,3,(((pl_indx-1)*3)+2))
        if max(x_cor(91:201,1))==0
            text(0,0,'No Cofiring');
        else
            bar(-10:100,x_cor(91:201,1))
            hold on
            line(-10:100,x_cor(91:201,2),'Color',[1 0 0]);
            line(-10:100,x_cor(91:201,4),'Color',[0.5 0.5 0.5]);
            line(-10:100,x_cor(91:201,5),'Color',[0.5 0.5 0.5]);
            line([0 0],[0 1.1*(max(max(x_cor(91:201,1:5))))],'Color',[0 1 0]);
            xlabel ('time after Pre spike (ms)');ylabel(' FR rate (post; Hz)');
            axis([-10, 100, 0, 1.1*(max(max(x_cor(91:201,1:5))))]);
            hold off
            %expand the graph
            pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi); clear pozi
            
            %percentages (red) and spikes (blue)
            subplot (5,3,(((pl_indx-1)*3)+3))
            hold on
            line(-10:100,x_cor(91:201,7),'Color',[1 0 0]);
            line(-10:100,x_cor(91:201,8),'Color',[0 0 1]);
            line([-10 100],[30 30],'Color',[0.5 0.5 0.5]);
            line([-10 100],[-30 -30],'Color',[0.5 0.5 0.5]);
            xlabel ('time after Pre spike (ms)');ylabel(' excess/missing %/# spikes');
            axis([-10, 100, (-20+(min(min(x_cor(91:201,7:8))))), (20+(max(max(x_cor(91:201,7:8)))))]);
            hold off
            %plot text results
            if char_intact{pr_indx,pt_indx,2}==1
                text(30,0,strcat('Cofiring:',num2str(char_intactQ{pr_indx,pt_indx,1})));
            end
            if char_intact{pr_indx,pt_indx,3}==1
                text(30,40,strcat('Excitation:',num2str(char_intactQ{pr_indx,pt_indx,2})));
            end
            if char_intact{pr_indx,pt_indx,4}==1
                text(30,-45,strcat('Inhibition:',num2str(char_intactQ{pr_indx,pt_indx,3})));
            end
            %expand the graph
            pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi); clear pozi
        end
        %increase the plot number
        if pl_indx<5
            pl_indx=pl_indx+1;
        else
            fg_name=strcat('fig',num2str(bs_num),num2str(bs_exp),'_CHARglass_intTHTonlyP',num2str(fg_indx));
            saveas(gcf,fg_name)
            clear fg_name
            pl_indx=1; fg_indx=fg_indx+1;
        end
        
    end; clear pt_indx
    fg_name=strcat('fig',num2str(bs_num),num2str(bs_exp),'_CHARglass_intTHTonlyP',num2str(fg_indx));
    saveas(gcf,fg_name)
    clear fg_name
end; clear pr_indx pl_indx fg_indx

datum=date;

%save the relevant data
writename=strcat(bs_name, bs_num, bs_exp, bs_dayname, '_unitcharGLASS_forintTHTonly.mat');
save(writename, 'datum', 'un_shnk','ch_uns', 'un_nums','char_GLunits','char_intact','char_xcorBACK','char_xcorFORW','t_rec','spks','spks_th','glassp_all','glassp_tht','-v7.3');


