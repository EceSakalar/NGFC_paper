%This program reads WAW.mat files and takes waveform average of the spikes of the juxta unit, analyses it
%and analyses its interaction with other cells recorded in the tetrodes 

clear all; %close all

%unit basic data entry
bs_name='B'; bs_num='181'; bs_exp='g'; bs_typ='sx'; bs_dayname='d07';

%give the silicon probe sampling rate (in Hz)

%channels in the spike2 files and originating tetrodes
ch_uns = [ ];
un_nums=[ ];                 %unit numbers (in kwik file)
ch_tet =[ ];

%glass unit data 
gl_uns = [20];      %channel to read
gl_nums= [1];       %unit number
gl_chan= [17];      %rec voltage channel
ch_rate = 20000;

%time limits of registration for glass units (periods start and end times)
%for each unit (column) and file (row)
gl_limS = [363];        %starts
gl_limE = [657];        %ends

%files to be included
fl_ord=[1];

%window for the spike extraction (in ms, one side), and resample rate for spike
wav_win = 3;
wav_rate = 100000;

tic
%unit stepper
for gl_indx=1:length(gl_uns);
    
    %file stepper
    for fl_indx=1:length(fl_ord)
        
        %read the mat file for spike trains
        fl_name=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'WAW.mat');        %file name
        ch_name=strcat(bs_name, bs_typ, bs_num, bs_exp,'_Ch', num2str(gl_uns(gl_indx)));
        ch_load=load(num2str(fl_name),num2str(ch_name));
        fl_spks=(ch_load.(ch_name).('times')); clear ch_load ch_name;
        fl_spks(fl_spks<=gl_limS(fl_indx,gl_indx))=[];
        fl_spks(fl_spks>=gl_limE(fl_indx,gl_indx))=[];
        
        %read the mat file for spike voltage recording
        ch_name=strcat(bs_name, bs_typ, bs_num, bs_exp,'_Ch', num2str(gl_chan(gl_indx)));
        ch_load=load(num2str(fl_name),num2str(ch_name));
        fl_chan=(ch_load.(ch_name).('values')); clear ch_load ch_name fl_name;
        
        %collect & store all spike times for the glass units
        if fl_indx==1;
            Gspks{gl_indx} = fl_spks-(gl_limS(fl_indx,gl_indx));
            t_incr = (gl_limE(fl_indx,gl_indx)-gl_limS(fl_indx,gl_indx));
        else
            Gspks{gl_indx} = [Gspks{gl_indx}; ((fl_spks-(gl_limS(fl_indx,gl_indx)))+t_incr)];
            t_incr = t_incr+((gl_limE(fl_indx,gl_indx)-gl_limS(fl_indx,gl_indx)));
        end
        
        if fl_indx == length(fl_ord)
           t_rec{gl_indx}=t_incr;
        end
        
        %collect and store spike shapes
        sp_file=zeros(length(fl_spks),((wav_win*(ch_rate/1000))*2)+1);
        if isempty(fl_spks)
            sp_file=[];
        else
            for sp_indx = 1:length(fl_spks)
                %extracting the wave on the original spike detected with 0.5 ms window
                spk_align = fl_chan((round(fl_spks(sp_indx)*ch_rate))-(0.5*(ch_rate/1000)):(round(fl_spks(sp_indx)*ch_rate))+(0.5*(ch_rate/1000)));
                %realign if not accurate
                [~,spk_maxpos] = max(spk_align);
                spk_corr=spk_maxpos-((0.5*(ch_rate/1000))+1); clear spk_align spk_maxpos;
                %extract the new corrected spike
                if ((round(fl_spks(sp_indx)*ch_rate))-(wav_win*(ch_rate/1000))+spk_corr)>0&&((round(fl_spks(sp_indx)*ch_rate))+(wav_win*(ch_rate/1000))+spk_corr)<size(fl_chan,1)
                    spk_wave = fl_chan(((round(fl_spks(sp_indx)*ch_rate))-(wav_win*(ch_rate/1000))+spk_corr):((round(fl_spks(sp_indx)*ch_rate))+(wav_win*(ch_rate/1000))+spk_corr));
                    clear spk_corr
                    sp_file(sp_indx,:) = spk_wave; clear spk_wave
                end
            end
        end
        clear sp_indx fl_spks
        
        %collect spike data over the files
        if fl_indx==1
            sp_all=sp_file;
        else
            sp_all=[sp_all; sp_file];
        end
        clear sp_file fl_chan
    end
    clear fl_indx
    
    char_GLunits{gl_indx,1} = 'glass';
    char_GLunits{gl_indx,2} = gl_indx;
    char_GLunits{gl_indx,3} = sp_all;
    char_GLunits{gl_indx,4} = mean(sp_all);
    %resampling & baseline correcting the spike (mean)
    char_GLunits{gl_indx,5}(1,:) = [(-1*wav_win):(1/wav_rate)*1000:wav_win];
    char_GLunits{gl_indx,5}(2,:) = interp1((-1*wav_win):(1/ch_rate)*1000:wav_win,mean(sp_all),(-1*wav_win):(1/wav_rate)*1000:wav_win)-mean(sp_all(:,1));
    char_GLunits{gl_indx,5}(3,:) = interp1((-1*wav_win):(1/ch_rate)*1000:wav_win,mean(sp_all),(-1*wav_win):(1/wav_rate)*1000:wav_win,'spline');
    
    %smoothing by ten samples, taking until 40%
    inter1 = (smooth(mean(sp_all),10));
    inter2 = inter1([[1:10:((round(length(mean(sp_all))*0.03333))*10)+1] floor(length(mean(sp_all))*0.8333)  floor(length(mean(sp_all))*0.916666) end]);
    inter3 = [(-1*wav_win):(1/ch_rate)*1000:wav_win];
    inter4 = inter3([[1:10:((round(length(mean(sp_all))*0.03333))*10)+1] floor(length(mean(sp_all))*0.8333)  floor(length(mean(sp_all))*0.916666) end]);
    char_GLunits{gl_indx,5}(4,:) = interp1(inter4,inter2,(-1*wav_win):(1/wav_rate)*1000:wav_win,'spline');
    char_GLunits{gl_indx,5}(5,:) = char_GLunits{gl_indx,5}(3,:) -  char_GLunits{gl_indx,5}(4,:); 
    clear sp_all inter1 inter2 inter3 inter4
    
    %spike charistics calculus
    %find the peak and the trough
    [aPeak,sPeak] = max(char_GLunits{gl_indx,5}(5,:));
    [aTrgh,sTrgh] = min(char_GLunits{gl_indx,5}(5,:));
    
    %find 'decay' to 0.9 before the peak and after the trough of spike
    sFirst = find(char_GLunits{gl_indx,5}(5,1:sPeak)<=aPeak*0.1, 1, 'last' );
    sLast = find(char_GLunits{gl_indx,5}(5,sTrgh:end)>=aTrgh*0.1, 1, 'first' );
    t1 = ((sFirst-sPeak)/wav_rate)*1000;
    t2 = (sLast/wav_rate)*1000;
    sym =(abs(t1)-abs(t2))/(abs(t1)+abs(t2));
    
    char_GLunits{gl_indx,6} = [t1 t2 aPeak sym];
    clear t1 t2 sym sPeak aPeak sFirst sLast sTrgh aTrgh
    
end
clear gl_indx t_incr
'first calculus finished'
toc

%generate autocorrelograms
for gl_indx=1:length(gl_uns);
    %make autocorrelogram via automatrix (+/- 100 ms, 1ms bins)
    %make an accomodating matrix
    a_mat=NaN(length(Gspks{gl_indx}),201);
    %do it straight
    if length(Gspks{gl_indx})<5000
        seg_mat0 = [zeros(100,size(Gspks{gl_indx},1)); repmat(Gspks{gl_indx},1,size(Gspks{gl_indx},1))-repmat(Gspks{gl_indx}',size(Gspks{gl_indx},1),1); zeros(100,size(Gspks{gl_indx},1))];
        for ka=0:-1:-200
            a_mat(:,((ka*-1)+1))=diag(seg_mat0,ka);
        end; clear ka seg_mat0
        
        %do it chunked
    else
        seg_num=ceil(length(Gspks{gl_indx})/4000);
        for seg_indx=1:seg_num
            if seg_indx==1
                seg_spks=Gspks{gl_indx}((((seg_indx-1)*4000)+1):(seg_indx*4000));
                seg_mat1=repmat(seg_spks,1,size(seg_spks,1))-repmat(seg_spks',size(seg_spks,1),1);
                seg_mat2=[zeros(100,4000); seg_mat1; zeros(100,4000)];
                seg_mat3=NaN(4000,201);
                for ka=0:-1:-200
                    seg_mat3(:,((ka*-1)+1))=diag(seg_mat2,ka);
                end; clear ka seg_spks seg_mat1 seg_mat2
                a_mat(1:3900,:)=seg_mat3(1:3900,:);
            elseif seg_indx==seg_num
                seg_spks=Gspks{gl_indx}((((seg_indx-1)*4000)-99):end);
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
                seg_spks=Gspks{gl_indx}((((seg_indx-1)*4000)-99):(seg_indx*4000));
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
    a_cor = (sum(hist(a_mat,-0.101:0.001:0.101),2))./(size(Gspks{gl_indx},1)).*1000;
    a_cor(1)=[]; a_cor(end)=[];a_cor(101)=0;
    %for the second column take only values >half maximum
    a_corS = a_cor; a_corS(a_cor<(max(a_cor)/2))=NaN;
    char_GLunits{gl_indx,7} = [a_cor a_corS];
    %calculate the mean rate, n and weighted avergae of >mean(rate) time values
    char_GLunits{gl_indx,8} = [length(Gspks{gl_indx})/t_rec{gl_indx} length(Gspks{gl_indx}) (nansum(abs([-100:1:100]').*a_corS))/nansum(a_corS)];
    clear a_mat a_cor a_corS
end
clear gl_indx 
'autocorrelograms done'
toc

%------Plot unit spike and autocorrelogram results
gl_indx=0; 
for pt_indx=1:ceil(length(gl_uns)/3)
    figure
    for ka=1:3
        gl_indx=gl_indx+1;
        if gl_indx<=length(gl_uns)
            %unit identity
            subplot (3,3,(((ka-1)*3)+1))
            axis off
            text(0,7,strcat(bs_name, bs_num, bs_exp,'(',bs_dayname,')'));
            text(0,5,strcat('unit:glass',num2str(gl_nums(gl_indx)),'in ch',num2str(gl_uns(gl_indx))));
            axis ([0 1 0 8])
            %spike and spike parameters
            subplot (3,3,(((ka-1)*3)+2))
            plot (char_GLunits{gl_indx,5}(1,:),char_GLunits{gl_indx,5}(2:5,:))
            set (gca,'ColorOrder',[0.2 0.2 0.2; 0.2 0.2 0.2;0.2 0.2 0.2;1 1 1]);
            axis tight
            xlabel ('time (ms)');ylabel('amplitude');
            text (-2.9 , (min(min(char_GLunits{gl_indx,5}(2:5,:)))+abs((min(min(char_GLunits{gl_indx,5}(2:5,:)))*0.5))) , strcat('t1=',num2str(char_GLunits{gl_indx,6}(1))));
            text (-2.9 , (min(min(char_GLunits{gl_indx,5}(2:5,:)))+abs((min(min(char_GLunits{gl_indx,5}(2:5,:)))*0.3))) , strcat('t2=',num2str(char_GLunits{gl_indx,6}(2))));
            text (-2.9 , (min(min(char_GLunits{gl_indx,5}(2:5,:)))+abs((min(min(char_GLunits{gl_indx,5}(2:5,:)))*0.1))) , strcat('sym=',num2str(char_GLunits{gl_indx,6}(4))));
            %expand the graph
            pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi)
            clear pozi
             
            %autocorrelogram and parameters
            subplot (3,3,(((ka-1)*3)+3))
            bar(-100:1:100,char_GLunits{gl_indx,7});
            xlabel ('time (ms)');ylabel('rate (Hz)');
            axis([-50, 50, 0, max(char_GLunits{gl_indx,7}(:,1))]);
            text (-49 , (max(char_GLunits{gl_indx,7}(:,1))*0.95) , strcat('rate=',num2str(char_GLunits{gl_indx,8}(1))));
            text (-49 , (max(char_GLunits{gl_indx,7}(:,1))*0.75) , strcat('n=',num2str(char_GLunits{gl_indx,8}(2))));
            text (-49 , (max(char_GLunits{gl_indx,7}(:,1))*0.55) , strcat('hub=',num2str(char_GLunits{gl_indx,8}(3))));
            %expand the graph
            pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi)
            clear pozi
            
        end
    end; clear ka
end; clear pt_indx gl_indx 


%read all spike trains for the silicon probe intaracting partners
for un_indx=1:length(ch_uns);
    %file stepper
    for fl_indx=1:length(fl_ord)
        %read the mat file for spike trains
        fl_name=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'WAW.mat');        %file name
        ch_name=strcat(bs_name, bs_typ, bs_num, bs_exp,'_Ch', num2str(ch_uns(un_indx)));
        ch_load=load(num2str(fl_name),num2str(ch_name));
        fl_spks=(ch_load.(ch_name).('times')); clear ch_load ch_name fl_name;
        
        for gl_indx=1:length(gl_uns);
            flg_spks=fl_spks;
            flg_spks(flg_spks<=gl_limS(fl_indx,gl_indx))=[];
            flg_spks(flg_spks>=gl_limE(fl_indx,gl_indx))=[];
            
            %collect & store all spike times for the glass units
            if fl_indx==1;
                Uspks{un_indx,gl_indx} = flg_spks-(gl_limS(fl_indx,gl_indx));
                t_incr = (gl_limE(fl_indx,gl_indx)-gl_limS(fl_indx,gl_indx));
            else
                Uspks{un_indx,gl_indx} = [Uspks{un_indx,gl_indx}; ((flg_spks-(gl_limS(fl_indx,gl_indx)))+t_incr)];
                t_incr = t_incr+((gl_limE(fl_indx,gl_indx)-gl_limS(fl_indx,gl_indx)));
            end; clear flg_spks
        end; clear gl_indx
    end; clear fl_indx fl_spks
end; clear un_indx t_incr

%------------------make cross correlograms via cross matrices     

for gl_indx=1:length(gl_uns);       %glass unit loop
    
    %loop for the xcorr partner
    for pt_indx = 1:length(ch_uns)
        
        %CALCULATE CROSS CORRELOGRAMS
        %if glass unit partner is n<5000 do it straight
        if length(Gspks{gl_indx})<5000;
            preSPKnum = length(Gspks{gl_indx});
            pstSPKnum = length(Uspks{pt_indx,gl_indx});
            x_mat = repmat(Uspks{pt_indx,gl_indx},1,preSPKnum) - repmat(Gspks{gl_indx}',pstSPKnum,1);
            x_mat(abs(x_mat)>1)=NaN;
            [~,V_indx]=find(~isnan(x_mat));
            M_indx=find(~isnan(x_mat));
            x_vect=x_mat(M_indx);
            clear M_indx x_mat
        %if presyn partner is big segment it!
        else 
            seg_num=ceil(length(Gspks{gl_indx})/4000);
            preSPKnum = length(Gspks{gl_indx});
            pstSPKnum = length(Uspks{pt_indx,gl_indx});
            x_vect = NaN;V_indx = NaN;
            for seg_indx=1:seg_num
                if seg_indx==1;
                    seg_spks = Gspks{gl_indx}((((seg_indx-1)*4000)+1):(seg_indx*4000));
                    seg_mat = repmat(Uspks{pt_indx,gl_indx},1,size(seg_spks,1)) - repmat(seg_spks',size(Uspks{pt_indx,gl_indx},1),1);
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
                    
                elseif seg_indx==seg_num;
                    seg_spks=Gspks{gl_indx}((((seg_indx-1)*4000)+1):end);
                    seg_mat = repmat(Uspks{pt_indx,gl_indx},1,size(seg_spks,1)) - repmat(seg_spks',size(Uspks{pt_indx,gl_indx},1),1);
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
                    seg_spks=Gspks{gl_indx}((((seg_indx-1)*4000)+1):(seg_indx*4000));
                    seg_mat = repmat(Uspks{pt_indx,gl_indx},1,size(seg_spks,1)) - repmat(seg_spks',size(Uspks{pt_indx,gl_indx},1),1);
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
                end
            end; clear seg_indx seg_num
            x_vect(1)=[];
            V_indx(1)=[];
        end
        
        %filling the xcor                
        x_cor = ((hist(x_vect,-0.101:0.001:0.101))./(preSPKnum).*1000)';
        x_cor(1)=[]; x_cor(end)=[];                
        char_GxcorPre{gl_indx,pt_indx}=x_cor;                
        
        %filling the 'mirror' x_cor
        x_cor = ((hist((x_vect.*-1),-0.101:0.001:0.101))./(pstSPKnum).*1000)';
        x_cor(1)=[]; x_cor(end)=[];
        char_GxcorPst{pt_indx,gl_indx}=x_cor;
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
        %SECOND layer of INTERACTION matrix
            %1 if INHIBITION
           
        if mean(((char_GxcorPre{gl_indx,pt_indx}(1:end,1)).*preSPKnum)./1000)<=1 
            char_GintactPre{gl_indx,pt_indx,1} = NaN;
            char_GintactPst{pt_indx,gl_indx,1} = NaN;
            
        else
            char_GintactPre{gl_indx,pt_indx,1} = 0;
            char_GintactPst{pt_indx,gl_indx,1} = 0;
            %SHUFFLING
            %generate 100 shuffling shu_xcor variable
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

            %forward direction (glass presynaptic)--------------------------------------
            x_cor = char_GxcorPre{gl_indx,pt_indx};                         %xcor
            x_cor(:,2) = (mean(shu_xcor,1)./(preSPKnum)).*1000;             %expected
            x_cor(:,3) = (std(shu_xcor,1,1)./(preSPKnum)).*1000;            %STD expectation
            x_cor(:,4) = x_cor(:,2)+(3.3.*(x_cor(:,3)));                    %confidence up
            x_cor(:,5) = x_cor(:,2)-(3.3.*(x_cor(:,3)));                    %confidence low
            x_cor(:,6) = NaN(201,1);                                        %over or under confidence
            x_cor(x_cor(:,1)>x_cor(:,4),6)=1;
            x_cor(x_cor(:,1)<x_cor(:,5),6)=-1;
            x_cor(:,7) = (((x_cor(:,1))./(x_cor(:,2)))-1)*100;              %excess or missing (percent of expected)
            x_cor(:,8) = (((x_cor(:,1))-(x_cor(:,2))).*preSPKnum)./1000;    %excess or missing (spikes)
            char_GxcorPre{gl_indx,pt_indx}=x_cor;
            %char_GxcorjitPre{gl_indx,pt_indx}=shu_xcor;
            
            %check if any significant bins in the forward direction(for pre(glass) - post interaction)
            if nansum(abs(x_cor(101:end,6)))>0
                char_GintactPre{gl_indx,pt_indx,1} = 1;
                char_GintactPre{gl_indx,pt_indx,2} = 0;
                char_GintactPre{gl_indx,pt_indx,3} = 0;
                char_GintactPre{gl_indx,pt_indx,4} = 0;
                %check cofiring
                if x_cor(101,1)>x_cor(101,4) && x_cor(101,7)>30 && x_cor(101,8)>30;
                    char_GintactPre{gl_indx,pt_indx,2} = 1;
                end
                x_eval=x_cor(102:105,:);
                %check excitation
                if ~isempty(x_eval(x_eval(:,6)==1,6))
                    if max(x_eval(x_eval(:,6)==1,7))>=30 && sum(x_eval(x_eval(:,6)==1,8))>=30;
                        char_GintactPre{gl_indx,pt_indx,3} = 1;
                    end
                end
                %check inhibition
                if ~isempty(x_eval(x_eval(:,6)==-1,6))
                    if min(x_eval(x_eval(:,6)==-1,7))<=-30 && sum(x_eval(x_eval(:,6)==-1,8))<=-30;
                        char_GintactPre{gl_indx,pt_indx,4} = 1;
                    end
                end
                clear x_eval
                strcat('pairs pre ', num2str(un_nums(pt_indx)), 'to post ', num2str(gl_nums(gl_indx)))
                toc
            end
            
            %reverse direction (glass postsynaptic)--------------------------------------
            x_cor = char_GxcorPst{pt_indx,gl_indx};
            shu_xcorB = zeros(size(shu_xcor));
            for ka=1:201; shu_xcorB(:,202-ka)=shu_xcor(:,ka); end; clear shu_xcor ka
            x_cor(:,2) = (mean(shu_xcorB,1)./(pstSPKnum)).*1000;             %expected
            x_cor(:,3) = (std(shu_xcorB,1,1)./(pstSPKnum)).*1000;            %STD expectation
            x_cor(:,4) = x_cor(:,2)+(3.3.*(x_cor(:,3)));                    %confidence up
            x_cor(:,5) = x_cor(:,2)-(3.3.*(x_cor(:,3)));                    %confidence low
            x_cor(:,6) = NaN(201,1);                                        %over or under confidence
            x_cor(x_cor(:,1)>x_cor(:,4),6)=1;
            x_cor(x_cor(:,1)<x_cor(:,5),6)=-1;
            x_cor(:,7) = (((x_cor(:,1))./(x_cor(:,2)))-1)*100;              %excess or missing (percent of expected)
            x_cor(:,8) = (((x_cor(:,1))-(x_cor(:,2))).*pstSPKnum)./1000;    %excess or missing (spikes)
            char_GxcorPst{pt_indx,gl_indx}=x_cor;
            %char_GxcorjitPst{pt_indx,gl_indx}=shu_xcorB;        
            
            %check if any significant bins in the reverse direction (for pre(silicon) - post interaction)
            if nansum(abs(x_cor(101:end,6)))>0
                char_GintactPst{pt_indx,gl_indx,1} = 1;
                char_GintactPst{pt_indx,gl_indx,2} = 0;
                char_GintactPst{pt_indx,gl_indx,3} = 0;
                char_GintactPst{pt_indx,gl_indx,4} = 0;
                %check cofiring
                if x_cor(101,1)>x_cor(101,4) && x_cor(101,7)>30 && x_cor(101,8)>30;
                    char_GintactPst{pt_indx,gl_indx,2} = 1;
                end
                x_eval=x_cor(102:105,:);
                %check excitation
                if ~isempty(x_eval(x_eval(:,6)==1,6))
                    if max(x_eval(x_eval(:,6)==1,7))>=30 && sum(x_eval(x_eval(:,6)==1,8))>=30;
                        char_GintactPst{pt_indx,gl_indx,3} = 1;
                    end
                end
                %check inhibition
                if ~isempty(x_eval(x_eval(:,6)==-1,6))
                    if min(x_eval(x_eval(:,6)==-1,7))<=-30 && sum(x_eval(x_eval(:,6)==-1,8))<=-30;
                        char_GintactPst{pt_indx,gl_indx,4} = 1;
                    end
                end
                clear x_eval
                strcat('pairs pre ', num2str(gl_nums(gl_indx)), 'to post ', num2str(un_nums(pt_indx)))
                toc
            end;
            clear x_cor shu_xcorB 
        end
    end; clear pt_indx
end
clear gl_indx preSPKnum pstSPKnum 

%------------PLOT CROSSCORELOGRAMS

%plot index
fg_indx=1; pl_indx=1;

for gl_indx=1:length(gl_uns);       %glass unit loop
    %loop for the xcorr partner
    for pt_indx = 1:length(ch_uns)
        %Is the interaction to be plotted? - FORWARD
        if char_GintactPre{gl_indx,pt_indx,1}==1;
            %new figure?
            if pl_indx==1;
                figure
            end
            
            %units involved identity, interactions
            subplot (5,3,(((pl_indx-1)*3)+1))
            axis off
            text(0,11,strcat(bs_name, bs_num, bs_exp,'(',bs_dayname,')'));
            text(0,9,strcat('Pre unit: glass',num2str(gl_nums(gl_indx))));
            text(0,7,strcat('Post unit: shank',num2str(ch_tet(pt_indx)),' unit', num2str(un_nums(pt_indx))));
            text(0,5,strcat('Cofiring:',num2str(char_GintactPre{gl_indx,pt_indx,2})));
            text(0,3,strcat('Excitation:',num2str(char_GintactPre{gl_indx,pt_indx,3})));
            text(0,1,strcat('Inhibition:',num2str(char_GintactPre{gl_indx,pt_indx,4})));
            axis ([0 1 0 12])
            
            %crosscorelations, expected rates and confidence intervals from shuffling
            subplot (5,3,(((pl_indx-1)*3)+2))
            bar(-10:100,char_GxcorPre{gl_indx,pt_indx}(91:201,1))
            hold on
            line(-10:100,char_GxcorPre{gl_indx,pt_indx}(91:201,2),'Color',[1 0 0]);
            line(-10:100,char_GxcorPre{gl_indx,pt_indx}(91:201,4),'Color',[0.5 0.5 0.5]);
            line(-10:100,char_GxcorPre{gl_indx,pt_indx}(91:201,5),'Color',[0.5 0.5 0.5]);
            line([0 0],[0 1.1*(max(max(char_GxcorPre{gl_indx,pt_indx}(91:201,1:5))))],'Color',[0 1 0]);
            xlabel ('time after Pre spike (ms)');ylabel(' FR rate (Hz)');
            axis([-10, 100, 0, 1.1*(max(max(char_GxcorPre{gl_indx,pt_indx}(91:201,1:5))))]);
            hold off
            %expand the graph
            pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi); clear pozi
            
            %percentages (red) and spikes (blue)
            subplot (5,3,(((pl_indx-1)*3)+3))
            hold on
            line(-10:100,char_GxcorPre{gl_indx,pt_indx}(91:201,7),'Color',[1 0 0]);
            line(-10:100,char_GxcorPre{gl_indx,pt_indx}(91:201,8),'Color',[0 0 1]);
            line([-10 100],[30 30],'Color',[0.5 0.5 0.5]);
            line([-10 100],[-30 -30],'Color',[0.5 0.5 0.5]);
            xlabel ('time after Pre spike (ms)');ylabel(' excess/missing %/# spikes');
            axis([-10, 100, (-20+(min(min(char_GxcorPre{gl_indx,pt_indx}(91:201,7:8))))), (20+(max(max(char_GxcorPre{gl_indx,pt_indx}(91:201,7:8)))))]);
            hold off
            %expand the graph
            pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi); clear pozi
            
            %increase the plot number
            if pl_indx<5
                pl_indx=pl_indx+1;
            else
                pl_indx=1; fg_indx=fg_indx+1;
            end
        end
        
        %Is the interaction to be plotted? - REVERSE
        if char_GintactPst{pt_indx,gl_indx,1}==1;
            %new figure
            if pl_indx==1;
                figure
            end
            
            %units involved identity, interactions
            subplot (5,3,(((pl_indx-1)*3)+1))
            axis off
            text(0,11,strcat(bs_name, bs_num, bs_exp,'(',bs_dayname,')'));
            text(0,7,strcat('Post unit: glass',num2str(gl_nums(gl_indx))));
            text(0,9,strcat('Pre unit: shank',num2str(ch_tet(pt_indx)),' unit', num2str(un_nums(pt_indx))));
            text(0,5,strcat('Cofiring:',num2str(char_GintactPst{pt_indx,gl_indx,2})));
            text(0,3,strcat('Excitation:',num2str(char_GintactPst{pt_indx,gl_indx,3})));
            text(0,1,strcat('Inhibition:',num2str(char_GintactPst{pt_indx,gl_indx,4})));
            axis ([0 1 0 12])
            
            %crosscorelations, expected rates and confidence intervals from shuffling
            subplot (5,3,(((pl_indx-1)*3)+2))
            bar(-10:100,char_GxcorPst{pt_indx,gl_indx}(91:201,1))
            hold on
            line(-10:100,char_GxcorPst{pt_indx,gl_indx}(91:201,2),'Color',[1 0 0]);
            line(-10:100,char_GxcorPst{pt_indx,gl_indx}(91:201,4),'Color',[0.5 0.5 0.5]);
            line(-10:100,char_GxcorPst{pt_indx,gl_indx}(91:201,5),'Color',[0.5 0.5 0.5]);
            line([0 0],[0 1.1*(max(max(char_GxcorPst{pt_indx,gl_indx}(91:201,1:5))))],'Color',[0 1 0]);
            xlabel ('time after Pre spike (ms)');ylabel(' FR rate (Hz)');
            axis([-10, 100, 0, 1.1*(max(max(char_GxcorPst{pt_indx,gl_indx}(91:201,1:5))))]);
            hold off
            %expand the graph
            pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi); clear pozi
            
            %percentages (red) and spikes (blue)
            subplot (5,3,(((pl_indx-1)*3)+3))
            hold on
            line(-10:100,char_GxcorPst{pt_indx,gl_indx}(91:201,7),'Color',[1 0 0]);
            line(-10:100,char_GxcorPst{pt_indx,gl_indx}(91:201,8),'Color',[0 0 1]);
            line([-10 100],[30 30],'Color',[0.5 0.5 0.5]);
            line([-10 100],[-30 -30],'Color',[0.5 0.5 0.5]);
            xlabel ('time after Pre spike (ms)');ylabel(' excess/missing %/# spikes');
            axis([-10, 100, (-20+(min(min(char_GxcorPst{pt_indx,gl_indx}(91:201,7:8))))), (20+(max(max(char_GxcorPst{pt_indx,gl_indx}(91:201,7:8)))))]);
            hold off
            %expand the graph
            pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi); clear pozi
            
            %increase the plot number
            if pl_indx<5
                pl_indx=pl_indx+1;
            else
                pl_indx=1; fg_indx=fg_indx+1;
            end
        end
    end; clear pt_indx
end; clear gl_indx pl_indx fg_indx
        
datum=date;

%save the relevant data
writename=strcat(bs_name, bs_num, bs_exp, bs_dayname, '_unitcharGLASS.mat');
save(writename, 'datum', 'ch_tet','ch_uns', 'un_nums','char_GintactPre','char_GintactPst','char_GLunits','char_GxcorPre','char_GxcorPst','t_rec','fl_ord','gl_chan','gl_nums','gl_uns','gl_chan','gl_limE','gl_limS','Gspks','Uspks');


