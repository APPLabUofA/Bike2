%SPECTRA PREPROCESSING AGAIN %do same for ERP if spectra looks improved, etc. 
parts = {'004' '005' '006' '007' '008' '009' '010' '011' '012' '013' '014' '016'};
%parts = {'004'};	
conditions = {'In' 'Out'};

filepath = ['M:\Data\bike\BikeOut\'];
exp = 'BikeOut';

eeg_thresh_1 = [-1000,1000];
eeg_thresh_2 = [-500,500];

if ~exist([filepath 'segments\']) %when preprocessing files, puts everything ina folder called segments, if not there, makes it
    mkdir([filepath 'segments_fft_JK\']);
end

high_pass = 0.1;
low_pass = 30;

epoch_size = [-1  2];
baseline = [-200    0];
eeg_thresh_1 = [-1000,1000];
eeg_thresh_2 = [-500,500];

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% i_part = 10;
% i_cond = 1;

%%
for i_part = 1:length(parts)
    
    for i_cond = 1:length(conditions)
        
%              
                 events = {'2'};
                 selection_cards = {'2'};
        %%%load data%%%
        disp(['Processing data for participant ' parts{i_part}]);
        disp(['Loading file: ' parts{i_part} '_' exp '_' conditions{i_cond} '.vhdr']);
        filename = [parts{i_part} '_' exp '_' conditions{i_cond} '.vhdr'];
        EEG = pop_loadbv(filepath, filename, [], []); % loading in EEG data
        savename = [parts{i_part} '_' exp '_' conditions{i_cond}];
        
        % get electrode locations
        EEG=pop_chanedit(EEG, 'load',{'M:\Analysis\Electrodelocs\Vamp_EOG_electrode_locs.ced' 'filetype' 'autodetect'});
        
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
        
        EEG = pop_epoch( EEG, events, epoch_size, 'newname',  sprintf('%s epochs' , exp), 'epochinfo', 'yes'); %Changed from [-.2 1] to [-1 2]. DR
        %S1 is standard, S2 = target....
        EEG = pop_rmbase( EEG, baseline);
        
        %    Artifact rejection, trials with range >500 uV
        EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)],eeg_thresh_1(1),eeg_thresh_1(2),EEG.xmin,EEG.xmax,0,1);
        
        %   EMCP occular correction
        temp_ocular = EEG.data(end-1:end,:,:); %to save the EYE data for after
        EEG = gratton_emcp(EEG,selection_cards,{'VEOG'},{'HEOG'}); %this assumes the eye channels are called this
        EEG.emcp.table %this prints out the regression coefficients
        EEG.data(end-1:end,:,:) = temp_ocular; %replace the eye data
        
        %    Artifact rejection, trials with range >250 uV
        EEG = pop_rmbase( EEG, baseline); %baseline again since this changed it
        EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)-2],eeg_thresh_2(1),eeg_thresh_2(2),EEG.xmin,EEG.xmax,0,1);
        
        EEG_Copy = EEG;
        
        EEG = pop_selectevent(EEG_Copy, 'type',str2num(events{1}),'renametype','Target','deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = pop_editset(EEG, 'setname',[savename '_fft_Target']);
        EEG = pop_saveset(EEG, 'filename',[savename '_fft_Target'],'filepath',[filepath 'segments_fft_JK\']);
        
    end
end