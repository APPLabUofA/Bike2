clear all
close all

exp = 'bike2';
subs = {'100' '101' '102'  '104' '106'  '108' '110'... 
        '114' '115' '116' '117' '118' '120' '121'...
        '122'  '126' '127' '129' '130' '131' '132' '133'...
         '135' '136'};

%subs = {'109'}; %to test on just one sub

nsubs = length(subs);
conds = {'sask'; '110st'; '83ave'};
conds_lab = {'Sask Drive'; '110 Street'; '83 Avenue'};
nconds = length(conds);
Pathname = 'M:\Data\Bike_lanes\';

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

%%
for i_sub = 1:nsubs
    for i_cond = 1:nconds
        
        Filename = [subs{i_sub} '_' exp '_' conds{i_cond}];
        EEG = pop_loadset('filename',[Filename '_Targets.set'],'filepath','M:\Data\Bike_lanes\segments_JK\');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        EEG = pop_loadset('filename',[Filename '_Standards.set'],'filepath','M:\Data\Bike_lanes\segments_JK');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
    end
end
eeglab redraw
%%
%subject erps
electrode = 13;%this is PZ in this electrode map
%electrode = [13,14]; %to average for electrodes 
erp_out = [];
for i_sub = 1:nsubs
    for i_cond = 1:nconds
        %average over trials (3rd dimension)
        erp_out(1,i_cond,i_sub) = length(ALLEEG(1+ 2*((i_sub-1)*nconds+(i_cond-1))).data(1,1,:)); %Targets
        erp_out(2,i_cond,i_sub) = length(ALLEEG(2+ 2*((i_sub-1)*nconds+(i_cond-1))).data(1,1,:)); %standards
    end
end

%%%get mean trial count for each condition%%%
%%%targets%%%
[mean(erp_out(1,1,:),3);mean(erp_out(1,2,:),3);mean(erp_out(1,3,:),3)]
t_sask = squeeze (erp_out(1,1,:));
t_110 = squeeze (erp_out(1,2,:));
t_83 = squeeze (erp_out(1,3,:));
%%%standards%%%
[mean(erp_out(2,1,:),3);mean(erp_out(2,2,:),3);mean(erp_out(2,3,:),3)]
s_sask = squeeze (erp_out(2,1,:));
s_110 = squeeze (erp_out(2,2,:));
s_83 = squeeze (erp_out(2,3,:));

%%%get SD trial count for each condition%%%
%%%targets%%%
[std(erp_out(1,1,:),[],3);std(erp_out(1,2,:),[],3);std(erp_out(1,3,:),[],3)]

%%%standards%%%
[std(erp_out(2,1,:),[],3);std(erp_out(2,2,:),[],3);std(erp_out(2,3,:),[],3)]


%%%get max and max trial count for each condition%%%
%%%targets%%%
[min(erp_out(1,1,:)),max(erp_out(1,1,:));min(erp_out(1,2,:)),max(erp_out(1,2,:));min(erp_out(1,3,:)),max(erp_out(1,3,:))]

%%%standards%%%
[min(erp_out(2,1,:)),max(erp_out(2,1,:));min(erp_out(2,2,:)),max(erp_out(2,2,:));min(erp_out(2,3,:)),max(erp_out(2,3,:))]