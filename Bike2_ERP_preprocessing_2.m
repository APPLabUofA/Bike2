%ERP PREPROCESSING AGAIN
ccc
parts = {'100' '101' '102' '103' '104' '106' '107' '108' '110'... 
        '113' '114' '115' '116' '117' '118' '119' '120' '121'...
        '122' '123' '126' '127' '129' '130' '131' '132' '133'...
        '134' '135' '136'};
% parts = {'125' '126'	'127'	'129' '131'	'132'};

 %parts = {'101' '102'};	%testing from the two trigger types

conditions = {'sask' '110st' '83ave'};

filepath = ['M:\Data\Bike_lanes\'];
exp = 'Bike2';

eeg_thresh_1 = [-1000,1000];
eeg_thresh_2 = [-500,500];

if ~exist([filepath 'segments\']) %when preprocessing files, puts everything ina folder called segments, if not there, makes it
    mkdir([filepath 'segments_JK\']);
end

high_pass = 0.1;
low_pass = 30;

epoch_size = [-.2  1];
baseline = [-200    0];
eeg_thresh_1 = [-1000,1000];
eeg_thresh_2 = [-500,500];

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


%%
for i_part = 1:length(parts)
    
    for i_cond = 1:length(conditions)
               
                 events = {'1' '2'};
                 selection_cards = {'1' '2'};
                 
        %%%load data%%%
        disp(['Processing data for participant ' parts{i_part}]);
        disp(['Loading file: ' parts{i_part} '_' exp '_' conditions{i_cond} '.vhdr']);
        filename = [parts{i_part} '_' exp '_' conditions{i_cond} '.vhdr'];
        EEG = pop_loadbv(filepath, filename, [], []); % loading in EEG data
        savename = [parts{i_part} '_' exp '_' conditions{i_cond}];
        
        % get electrode locations
        EEG=pop_chanedit(EEG, 'load',{'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced' 'filetype' 'autodetect'});
        
        % arithmetically rereference to linked mastoid
        for x=1:EEG.nbchan-2
            EEG.data(x,:) = (EEG.data(x,:)-((EEG.data(EEG.nbchan-2,:))*.5));
        end
        
                %Filter the data with low pass of 30
                EEG = pop_eegfilt( EEG, high_pass, 0, [], 0);  %high pass filter
                EEG = pop_eegfilt( EEG, 0, low_pass, [], 0);  %low pass filter
        
        
        all_events = length(EEG.event);
        
        
        for i_event = 2:all_events %%% Why start from 2
            if strcmp(EEG.event(i_event).type, 'boundary')
                continue
            else
                EEG.event(i_event).type = num2str(str2num((EEG.event(i_event).type(2:end))));
            end
        end
        
        %epoch
        
        EEG = pop_epoch( EEG, events, epoch_size, 'newname',  sprintf('%s epochs' , exp), 'epochinfo', 'yes'); 
        %S1 is standard, S2 = target....
        EEG = pop_rmbase( EEG, baseline);
        
        %    Artifact rejection, trials with range >1000 uV
        EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)],eeg_thresh_1(1),eeg_thresh_1(2),EEG.xmin,EEG.xmax,0,1);
        
        %   EMCP occular correction
        temp_ocular = EEG.data(end-1:end,:,:); %to save the EYE data for after
        EEG = gratton_emcp(EEG,selection_cards,{'VEOG'},{'HEOG'}); %this assumes the eye channels are called this
        EEG.emcp.table %this prints out the regression coefficients
        EEG.data(end-1:end,:,:) = temp_ocular; %replace the eye data
        
        %    Artifact rejection, trials with range >500 uV
        EEG = pop_rmbase( EEG, baseline); %baseline again since this changed it
        EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)-2],eeg_thresh_2(1),eeg_thresh_2(2),EEG.xmin,EEG.xmax,0,1);
        
        EEG_Copy = EEG;
        
        EEG = pop_selectevent(EEG_Copy, 'type',str2num(events{2}),'renametype','Targets','deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = pop_editset(EEG, 'setname',[savename '_Targets']);
        EEG = pop_saveset(EEG, 'filename',[savename '_Targets'],'filepath',[filepath 'segments_JK\']);
        
        EEG = pop_selectevent(EEG_Copy, 'type',str2num(events{1}),'renametype','Standards','deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = pop_editset(EEG, 'setname',[savename '_Standards']);
        EEG = pop_saveset(EEG, 'filename',[savename '_Standards'],'filepath',[filepath 'segments_JK\']);
        
        
        
    end
end