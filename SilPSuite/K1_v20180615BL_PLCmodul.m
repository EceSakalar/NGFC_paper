%This program analyses the 2D place-dependent firing in all the cells

clear all; close all

bs_name='B'; bs_num='178'; bs_exp='a'; bs_typ='sx';                 %unit basic data entry

ch_n=16;                                                            %number of channels
ch_rate=2000;                                                       %give the silicon probe sampling rate (in Hz)

un_chans=[19:28];                                                      %channels containing the spike times for the units
un_nums=[5 35 37 43 46 47 54 2 11 16];                                                        %unit numbers (in kwik file)
un_shank =[repmat([1],1,7) repmat([2],1,3)];                       %shanks from which unit was isolated

ch_def=[];                                                          %defective channel (one is standing for the ventralmost, and the linear number of channel is required not the identifier
fl_ord=[0];                                                         %file extension numbers belonging to the present cell
