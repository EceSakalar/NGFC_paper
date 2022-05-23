%This program reads original .dat files and takes waveform averages 
%for the optimal channels (where the maximum amplitude is observed, calculates auto and cross corelograms

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
        acc_chan = acc_data(:,(ch_max(un_indx))); clear acc_data;
        
        sp_file=zeros(length(acc_spks),((wav_win*(ch_rate/1000))*2)+1);
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
                    acc_wave = acc_chan(((round(acc_spks(sp_indx)*ch_rate))-(wav_win*(ch_rate/1000))+acc_corr):((round(acc_spks(sp_indx)*ch_rate))+(wav_win*(ch_rate/1000))+acc_corr));
                    clear acc_corr
                    sp_file(sp_indx,:) = acc_wave; clear acc_wave
                end
            end
        end
        clear sp_indx acc_spks
        
        %collect spike data over the files
        if fl_indx==1
            sp_all=sp_file;
        else
            sp_all=[sp_all; sp_file];
        end
        clear sp_file acc_chan
    end
    clear fl_indx t_incr
    
    char_units{un_indx,1} = ch_tet(un_indx);
    char_units{un_indx,2} = un_indx;
    char_units{un_indx,3} = sp_all;
    char_units{un_indx,4} = mean(sp_all);
    
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
    
    %spike charistics calculus
    [aPeak,sPeak] = min(char_units{un_indx,5}(5,:));
    %find 'decay' to 0.9 before and after the spike trough
    sFirst = find(char_units{un_indx,5}(5,1:sPeak)>=aPeak*0.1, 1, 'last' );
    sLast = find(char_units{un_indx,5}(5,sPeak:end)>=aPeak*0.1, 1, 'first' ); sLast=sLast+(wav_win*wav_rate/1000);
    t1 = ((sFirst-sPeak)/wav_rate)*1000;
    t2 = ((sLast-sPeak)/wav_rate)*1000;
    sym =(abs(t1)-abs(t2))/(abs(t1)+abs(t2));
    
    char_units{un_indx,6} = [t1 t2 aPeak sym];
    clear t1 t2 sym sPeak aPeak sFirst sLast 
    
end
clear un_indx
'first calculus finished'
toc

%generate autocorrelograms
for un_indx=1:length(ch_uns);
    %make autocorrelogram via automatrix (+/- 100 ms, 1ms bins)
    %make an accomodating matrix
    a_mat=NaN(length(spks{un_indx}),201);
    %do it straight
    if length(spks{un_indx})<5000
        seg_mat0 = [zeros(100,size(spks{un_indx},1)); repmat(spks{un_indx},1,size(spks{un_indx},1))-repmat(spks{un_indx}',size(spks{un_indx},1),1); zeros(100,size(spks{un_indx},1))];
        for ka=0:-1:-200
            a_mat(:,((ka*-1)+1))=diag(seg_mat0,ka);
        end; clear ka seg_mat0
        
        %do it chunked
    else
        seg_num=ceil(length(spks{un_indx})/4000);
        for seg_indx=1:seg_num
            if seg_indx==1
                seg_spks=spks{un_indx}((((seg_indx-1)*4000)+1):(seg_indx*4000));
                seg_mat1=repmat(seg_spks,1,size(seg_spks,1))-repmat(seg_spks',size(seg_spks,1),1);
                seg_mat2=[zeros(100,4000); seg_mat1; zeros(100,4000)];
                seg_mat3=NaN(4000,201);
                for ka=0:-1:-200
                    seg_mat3(:,((ka*-1)+1))=diag(seg_mat2,ka);
                end; clear ka seg_spks seg_mat1 seg_mat2
                a_mat(1:3900,:)=seg_mat3(1:3900,:);
            elseif seg_indx==seg_num
                seg_spks=spks{un_indx}((((seg_indx-1)*4000)-99):end);
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
                seg_spks=spks{un_indx}((((seg_indx-1)*4000)-99):(seg_indx*4000));
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
    a_cor = (sum(hist(a_mat,-0.101:0.001:0.101),2))./(size(spks{un_indx},1)).*1000;
    a_cor(1)=[]; a_cor(end)=[];a_cor(101)=0;
    %for the second column take only values > half of the maximum
    a_corS = a_cor; a_corS(a_cor<(max(a_cor)/2))=NaN;
    char_units{un_indx,7} = [a_cor a_corS];
    %calculate the mean rate, n and weighted avergae of >mean(rate) time values
    char_units{un_indx,8} = [length(spks{un_indx})/t_rec length(spks{un_indx}) (nansum(abs([-100:1:100]').*a_corS))/nansum(a_corS)];
    clear a_mat a_cor a_corS
end
'autocorrelograms done'
toc

%------Plot unit spike and autocorrelogram results
un_indx=0; 
for pt_indx=1:ceil(length(ch_uns)/5)
    
    figure
    for ka=1:5
        un_indx=un_indx+1;
        if un_indx<=length(ch_uns)
            %unit identity
            subplot (5,3,(((ka-1)*3)+1))
            axis off
            text(0,7,strcat(bs_name, bs_num, bs_exp,'(',bs_dayname,')'));
            text(0,5,strcat('unit:t',num2str(ch_tet(un_indx)),'u',num2str(un_nums(un_indx)),'in ch',num2str(ch_uns(un_indx))));
            axis ([0 1 0 8])
            %spike and spike parameters
            subplot (5,3,(((ka-1)*3)+2))
            plot (char_units{un_indx,5}(1,:),char_units{un_indx,5}(2:5,:))
            set (gca,'ColorOrder',[0.2 0.2 0.2; 0.2 0.2 0.2;0.2 0.2 0.2;1 1 1]);
            axis tight
            xlabel ('time (ms)');ylabel('amplitude');
            text (-2.9 , (min(min(char_units{un_indx,5}(2:5,:)))+abs((min(min(char_units{un_indx,5}(2:5,:)))*0.5))) , strcat('t1=',num2str(char_units{un_indx,6}(1))));
            text (-2.9 , (min(min(char_units{un_indx,5}(2:5,:)))+abs((min(min(char_units{un_indx,5}(2:5,:)))*0.3))) , strcat('t2=',num2str(char_units{un_indx,6}(2))));
            text (-2.9 , (min(min(char_units{un_indx,5}(2:5,:)))+abs((min(min(char_units{un_indx,5}(2:5,:)))*0.1))) , strcat('sym=',num2str(char_units{un_indx,6}(4))));
            %expand the graph
            pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi)
            clear pozi 
             
            %autocorrelogram and parameters
            subplot (5,3,(((ka-1)*3)+3))
            bar(-100:1:100,char_units{un_indx,7});
            xlabel ('time (ms)');ylabel('rate (Hz)');
            axis([-50, 50, 0, max(char_units{un_indx,7}(:,1))]);
            text (-49 , (max(char_units{un_indx,7}(:,1))*0.95) , strcat('rate=',num2str(char_units{un_indx,8}(1))));
            text (-49 , (max(char_units{un_indx,7}(:,1))*0.75) , strcat('n=',num2str(char_units{un_indx,8}(2))));
            text (-49 , (max(char_units{un_indx,7}(:,1))*0.55) , strcat('hub=',num2str(char_units{un_indx,8}(3))));
            %expand the graph
            pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi)
            clear pozi
            
        end
    end; clear ka
end; clear pt_indx un_indx 
    
%make cross correlograms via cross matrices     
for un_indx=1:length(ch_uns);
    char_intact{un_indx,length(ch_uns)+1}=strcat('un1:t',num2str(ch_tet(un_indx)),'u',num2str(un_nums(un_indx)));
    
    %loop for the xcorr partner
    for pt_indx = un_indx+1:length(ch_uns)
        
        %if presynaptic partner is n<5000 do it straight
        if length(spks{un_indx})<5000;
            preSPKnum = length(spks{un_indx});
            pstSPKnum = length(spks{pt_indx});
            x_mat = repmat(spks{pt_indx},1,size(spks{un_indx},1)) - repmat(spks{un_indx}',size(spks{pt_indx},1),1);
            x_mat(abs(x_mat)>1)=NaN;
            [~,V_indx]=find(~isnan(x_mat));
            M_indx=find(~isnan(x_mat));
            x_vect=x_mat(M_indx);
            clear M_indx x_mat
            
        %if presyn partner is big segment it!
        else 
            seg_num=ceil(length(spks{un_indx})/4000);
            preSPKnum = length(spks{un_indx});
            pstSPKnum = length(spks{pt_indx});
            x_vect = NaN;V_indx = NaN;
            for seg_indx=1:seg_num
                if seg_indx==1;
                    seg_spks = spks{un_indx}((((seg_indx-1)*4000)+1):(seg_indx*4000));
                    seg_mat = repmat(spks{pt_indx},1,size(seg_spks,1)) - repmat(seg_spks',size(spks{pt_indx},1),1);
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
                    seg_spks=spks{un_indx}((((seg_indx-1)*4000)+1):end);
                    seg_mat = repmat(spks{pt_indx},1,size(seg_spks,1)) - repmat(seg_spks',size(spks{pt_indx},1),1);
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
                    seg_spks=spks{un_indx}((((seg_indx-1)*4000)+1):(seg_indx*4000));
                    seg_mat = repmat(spks{pt_indx},1,size(seg_spks,1)) - repmat(seg_spks',size(spks{pt_indx},1),1);
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
        char_xcor{un_indx,pt_indx}=x_cor;                
        
        %filling the 'mirror' x_cor
        x_cor = ((hist((x_vect.*-1),-0.101:0.001:0.101))./(pstSPKnum).*1000)';
        x_cor(1)=[]; x_cor(end)=[];
        char_xcor{pt_indx,un_indx}=x_cor;
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
            
        if mean(((char_xcor{un_indx,pt_indx}(102:end,1)).*(char_units{un_indx,8}(1,2)))./1000)<=1 && mean(((char_xcor{un_indx,pt_indx}(1:100,1)).*(char_units{un_indx,8}(1,2)))./1000)<=1
            char_intact{un_indx,pt_indx,1} = NaN;
            char_intact{pt_indx,un_indx,1} = NaN;
            
        else
            char_intact{un_indx,pt_indx,1} = 0;
            char_intact{pt_indx,un_indx,1} = 0;
            %generate shuffling shu_xcor variable
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
            x_cor = char_xcor{un_indx,pt_indx};
            x_cor(:,2) = (mean(shu_xcor,1)./(size(spks{un_indx},1))).*1000;
            x_cor(:,3) = (std(shu_xcor,1,1)./(size(spks{un_indx},1))).*1000;
            x_cor(:,4) = x_cor(:,2)+(3.3.*(x_cor(:,3)));
            x_cor(:,5) = x_cor(:,2)-(3.3.*(x_cor(:,3)));
            x_cor(:,6) = NaN(201,1);                                        %over or under confidence
            x_cor(x_cor(:,1)>x_cor(:,4),6)=1;
            x_cor(x_cor(:,1)<x_cor(:,5),6)=-1;
            x_cor(:,7) = (((x_cor(:,1))./(x_cor(:,2)))-1)*100;              %excess or missing (percent of expected)
            x_cor(:,8) = (((x_cor(:,1))-(x_cor(:,2))).*preSPKnum)./1000;    %excess or missing (spikes)
            char_xcor{un_indx,pt_indx}=x_cor;
            %char_xcorjit{un_indx,pt_indx}=shu_xcor;
            
            %check if any significant bins in the forward direction(for pre - post interaction)
            %for significant interactions store the quantification in INTACTQ matrix layers 1 (cofir), 2 (excitation, and 3 (inhibition) 
            if nansum(abs(x_cor(101:end,6)))>0
                char_intact{un_indx,pt_indx,1} = 1;
                char_intact{un_indx,pt_indx,2} = 0;
                char_intact{un_indx,pt_indx,3} = 0;
                char_intact{un_indx,pt_indx,4} = 0;
                %check cofiring
                if x_cor(101,1)>x_cor(101,4) && x_cor(101,7)>30 && x_cor(101,8)>30;
                    char_intact{un_indx,pt_indx,2} = 1;
                    char_intactQ{un_indx,pt_indx,1} = (x_cor(101,8))/preSPKnum;
                end
                x_eval=x_cor(102:105,:);
                %check excitation
                if ~isempty(x_eval(x_eval(:,6)==1,6))
                    if max(x_eval(x_eval(:,6)==1,7))>=30 && sum(x_eval(x_eval(:,6)==1,8))>=30;
                         char_intact{un_indx,pt_indx,3} = 1;
                         char_intactQ{un_indx,pt_indx,2} = sum(x_eval(x_eval(:,6)==1,8))/preSPKnum;
                    end
                end
                %check inhibition
                if ~isempty(x_eval(x_eval(:,6)==-1,6))
                    if min(x_eval(x_eval(:,6)==-1,7))<=-30 && sum(x_eval(x_eval(:,6)==-1,8))<=-30;
                        char_intact{un_indx,pt_indx,4} = 1;
                        char_intactQ{un_indx,pt_indx,3} = sum(x_eval(x_eval(:,6)==-1,8))/preSPKnum;
                    end
                end
                clear x_eval
                strcat('pairs pre ', num2str(un_nums(pt_indx)), 'to post ', num2str(un_nums(un_indx)))
                toc
            end
            
            %inverted direction--------------------------------------
            x_cor = char_xcor{pt_indx,un_indx};
            temp1 = (mean(shu_xcor,1)./(size(spks{pt_indx},1))).*1000;
            temp2 = temp1;
            for ka=1:length(temp1)
                temp2(ka) = temp1(length(temp1)-ka+1);
            end
            x_cor(:,2) = temp2; clear temp1 temp2 ka
            
            temp1 = (std(shu_xcor,1,1)./(size(spks{pt_indx},1))).*1000;
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
            char_xcor{pt_indx,un_indx}=x_cor;
            
            %check if any significant bins in the reverse direction(for pre - post interaction)
            %for significant interactions store the quantification in INTACTQ matrix layers 1 (cofir), 2 (excitation, and 3 (inhibition) 
            if nansum(abs(x_cor(101:end,6)))>0
                char_intact{pt_indx,un_indx,1} = 1;
                char_intact{pt_indx,un_indx,2} = 0;
                char_intact{pt_indx,un_indx,3} = 0;
                char_intact{pt_indx,un_indx,4} = 0;
                %check cofiring
                if x_cor(101,1)>x_cor(101,4) && x_cor(101,7)>30 && x_cor(101,8)>30;
                    char_intact{pt_indx,un_indx,2} = 1;
                    char_intactQ{pt_indx,un_indx,1} = (x_cor(101,8))/pstSPKnum;
                end
                x_eval=x_cor(102:105,:);
                %check excitation
                if ~isempty(x_eval(x_eval(:,6)==1,6))
                    if max(x_eval(x_eval(:,6)==1,7))>=30 && sum(x_eval(x_eval(:,6)==1,8))>=30;
                         char_intact{pt_indx,un_indx,3} = 1;
                         char_intactQ{pt_indx,un_indx,2} = sum(x_eval(x_eval(:,6)==1,8))/pstSPKnum;
                    end
                end
                %check inhibition
                if ~isempty(x_eval(x_eval(:,6)==-1,6))
                    if min(x_eval(x_eval(:,6)==-1,7))<=-30 && sum(x_eval(x_eval(:,6)==-1,8))<=-30;
                        char_intact{pt_indx,un_indx,4} = 1;
                        char_intactQ{pt_indx,un_indx,3} = sum(x_eval(x_eval(:,6)==-1,8))/pstSPKnum;
                    end
                end
                clear x_eval
                strcat('pairs pre ', num2str(un_nums(un_indx)), 'to post ', num2str(un_nums(pt_indx)))
                toc
            end
            clear x_cor shu_xcor
        end
    end; clear pt_indx
end
clear un_indx preSPKnum pstSPKnum 


%------------PLOT CROSSCORELOGRAMS

%plot index
fg_indx=1; pl_indx=1;

for pr_indx=1:length(ch_uns);       %glass unit loop
    %loop for the xcorr partner
    for pt_indx = 1:length(ch_uns)
        
        %Is the interaction to be plotted? - FORWARD
        %if max spike count is below 15 it is not plotted 
        if pr_indx~=pt_indx && char_intact{pr_indx,pt_indx,1}==1 && max(((char_xcor{pr_indx,pt_indx}(91:201,1)).*length(spks{pr_indx}))./1000)>=15;
            
            %new figure?
            if pl_indx==1;
                figure
            end
            
            %units involved identity, interactions
            subplot (5,3,(((pl_indx-1)*3)+1))
            axis off
            text(0,11,strcat(bs_name, bs_num, bs_exp,'(',bs_dayname,')'));
            text(0,9,strcat('Pre unit: shank',num2str(ch_tet(pr_indx)),' unit', num2str(un_nums(pr_indx))));
            text(0,7,strcat('Post unit: shank',num2str(ch_tet(pt_indx)),' unit', num2str(un_nums(pt_indx))));
            text(0,5,strcat('Cofiring:',num2str(char_intact{pr_indx,pt_indx,2})));
            text(0,3,strcat('Excitation:',num2str(char_intact{pr_indx,pt_indx,3})));
            text(0,1,strcat('Inhibition:',num2str(char_intact{pr_indx,pt_indx,4})));
            axis ([0 1 0 12])
            
            %crosscorelations, expected rates and confidence intervals from shuffling
            subplot (5,3,(((pl_indx-1)*3)+2))
            bar(-10:100,char_xcor{pr_indx,pt_indx}(91:201,1))
            hold on
            line(-10:100,char_xcor{pr_indx,pt_indx}(91:201,2),'Color',[1 0 0]);
            line(-10:100,char_xcor{pr_indx,pt_indx}(91:201,4),'Color',[0.5 0.5 0.5]);
            line(-10:100,char_xcor{pr_indx,pt_indx}(91:201,5),'Color',[0.5 0.5 0.5]);
            line([0 0],[0 1.1*(max(max(char_xcor{pr_indx,pt_indx}(91:201,1:5))))],'Color',[0 1 0]);
            xlabel ('time after Pre spike (ms)');ylabel(' FR rate (post; Hz)');
            axis([-10, 100, 0, 1.1*(max(max(char_xcor{pr_indx,pt_indx}(91:201,1:5))))]);
            hold off
            %expand the graph
            pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi); clear pozi
            
            %percentages (red) and spikes (blue)
            subplot (5,3,(((pl_indx-1)*3)+3))
            hold on
            line(-10:100,char_xcor{pr_indx,pt_indx}(91:201,7),'Color',[1 0 0]);
            line(-10:100,char_xcor{pr_indx,pt_indx}(91:201,8),'Color',[0 0 1]);
            line([-10 100],[30 30],'Color',[0.5 0.5 0.5]);
            line([-10 100],[-30 -30],'Color',[0.5 0.5 0.5]);
            xlabel ('time after Pre spike (ms)');ylabel(' excess/missing %/# spikes');
            axis([-10, 100, (-20+(min(min(char_xcor{pr_indx,pt_indx}(91:201,7:8))))), (20+(max(max(char_xcor{pr_indx,pt_indx}(91:201,7:8)))))]);
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
            
            %increase the plot number
            if pl_indx<5
                pl_indx=pl_indx+1;
            else
                pl_indx=1; fg_indx=fg_indx+1;
            end
        end
    end; clear pt_indx
end; clear pr_indx pl_indx fg_indx

datum=date;

%save the relevant data
writename=strcat(bs_name, bs_num, bs_exp, bs_dayname, '_unitcharSILPR.mat');
save(writename, 'datum', 'ch_tet','ch_uns', 'un_nums','char_units','char_intact','char_intactQ','char_xcor','t_rec','fl_ord','spks');


