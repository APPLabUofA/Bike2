function BikeOut_ERSP_2_Plot(ALLEEG, conds, EEG, elec_names, electrode, electrode_loc,...
    ersp, exp, Filename, freqs, itc, Pathname, perms, powbase, subs, times,...
    trialevent)
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%                             INFORMATION
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% SPECTOGRAM
% A spectogram is a 3d figure that plots time on the x-axis, frequency on the 
% y-axis, and shows you the power or phase-locking value for each point. 
% We compute spectograms if we have power and phase information, averaged 
% across trials, for at least one electrode. 
% This can help us understand the changes of power and phase throughout the 
% trial.

% TOPOPLOT
% A Topoplot is a graph of the distribution of power or phase-locking magnitude
% across the scalp, averaging a range of both time and frequency.
% We make a topoplot if we have power or phase-locking information, averaged 
% across the trials, for every electrode.
% This can help us understand the scalp distribution of power or phase-locking 
% at a crucial time and at the frequency of maximal effect.

% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% Variables working with:
% ersp(i_sub,i_cond,i_perm,i_chan,:,:)
% itc(i_sub,i_cond,i_perm,i_chan,:,:)
% powbase,times,freqs

% The variables ersp and itc will be a 6D variable: 
% (participants x conditions x events x electrodes x frequencies x timepoints)


eeglab redraw

% /////////////////////////////////////////////////////////////////////////
%%     ERSP plots averaged over subjects and channels: Out vs In
% /////////////////////////////////////////////////////////////////////////

for i_event = 1:length(trialevent)
    
    % (participants x conditions x events x electrodes x frequencies x timepoints)
%    in_ersp = squeeze(mean(mean(ersp(:,1,i_event,:,:,:),4),1));
%    out_ersp = squeeze(mean(mean(ersp(:,2,i_event,:,:,:),4),1));
    pref_ersp = squeeze(mean(mean(ersp(:,1,i_event,:,:,:),4),1));
    non_pref_ersp = squeeze(mean(mean(ersp(:,2,i_event,:,:,:),4),1));


    figure; 
    CLim = [-1.5 1.5];
    colormap('jet')

    % Subplot 1: non-pref
    subplot(3,1,1); 
    imagesc(times,freqs,non_pref_ersp);
    title(['ERSP: Preferred; ' char(trialevent{i_event})]); 
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
    ylabel('Freq (Hz)');
    colorbar
    % Subplot 1: pref
    subplot(3,1,2); 
    imagesc(times,freqs,pref_ersp);
    title(['ERSP: Non-Preferred; ' char(trialevent{i_event})]); 
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
    ylabel('Freq (Hz)');
    colorbar
    % Subplot 1: out-in
    subplot(3,1,3); 
    imagesc(times,freqs,pref_ersp-non_pref_ersp); 
    title(['ERSP: Preferred-Non-preferred; ' char(trialevent{i_event})]);
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
    ylabel('Freq (Hz)'); xlabel('Time (ms)');
    colorbar

    % Overall subplot title
%     supertitle(['ERSP: ' char(trialevent{i_event})],'FontSize',12);

    clear non_pref_ersp pref_ersp

end

clear i_event

% /////////////////////////////////////////////////////////////////////////
%%       ITC plots averaged over subjects and channels: Out vs In
% /////////////////////////////////////////////////////////////////////////


for i_event = 1:length(trialevent)
    
    % The variable will be a 6D variable: 
    % (participants x conditions x events x electrodes x frequencies x timepoints)
    in_itc = squeeze(mean(mean(abs(itc(:,1,i_event,:,:,:)),4),1));
    out_itc = squeeze(mean(mean(abs(itc(:,2,i_event,:,:,:)),4),1));

    CLim = [-.5 .5];
    figure; 
    colormap('jet')

    % Subplot 1: out
    subplot(3,1,1); 
    imagesc(times,freqs,out_itc,CLim);
    title(['ITC: Out; ' char(trialevent{i_event})]); 
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
    ylabel('Freq (Hz)');
    colorbar
    % Subplot 1: in
    subplot(3,1,2); 
    imagesc(times,freqs,in_itc,CLim);
    title(['ITC: In; ' char(trialevent{i_event})]); 
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
    ylabel('Freq (Hz)');
    colorbar
    % Subplot 1: out-in
    subplot(3,1,3); 
    imagesc(times,freqs,out_itc-in_itc,CLim); 
    title(['ERSP: Out-In; ' char(trialevent{i_event})]);
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
    ylabel('Freq (Hz)'); xlabel('Time (ms)');
    colorbar

    %Overall subplot title
    % supertitle('ITC: Overall Average','FontSize',12)

    % close all
    % eeglab redraw
    clear out_itc in_itc  

end

clear i_event
    
% /////////////////////////////////////////////////////////////////////////
%%         ERSP plots averaged over subjects for each channel
% /////////////////////////////////////////////////////////////////////////


CLim = [-1.5 1.5];

for i_event = 1:length(trialevent)

    for ch_ersp = 1:length(electrode)

        % ERSP values by electrode
        % (participants x conditions x events x electrodes x frequencies x timepoints)
        pref_ersp_chan  =  squeeze(mean(ersp(:,1,i_event,ch_ersp,:,:),1));
        npref_ersp_chan = squeeze(mean(ersp(:,2,i_event,ch_ersp,:,:),1));

        figure; 
        colormap('jet')

        % Subplot 1: out
        subplot(3,1,1); 
        imagesc(times,freqs,npref_ersp_chan,CLim);
        title({['ERSP: ' char(trialevent{i_event}) ', NonPref']; char(elec_names(ch_ersp))}); 
        set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
        ylabel('Freq (Hz)');
        colorbar
        % Subplot 1: in
        subplot(3,1,2); 
        imagesc(times,freqs,pref_ersp_chan,CLim);
        title({['ERSP: ' char(trialevent{i_event}) ', Pref']; char(elec_names(ch_ersp))}); 
        set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
        ylabel('Freq (Hz)');
        colorbar
        % Subplot 1: out-in
        subplot(3,1,3); 
        imagesc(times,freqs,npref_ersp_chan-pref_ersp_chan,CLim); 
        title({['ERSP: ' char(trialevent{i_event}) ', np-p']; char(elec_names(ch_ersp))});
        set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
        ylabel('Freq (Hz)'); xlabel('Time (ms)');
        colorbar

        %Overall subplot title
    %     supertitle(['ERSP: ' char(elec_names(ch_ersp))],'FontSize',12)

        clear in_ersp_chan out_ersp_chan 
    end

    clear ch_ersp 

end

clear i_event

% close all
% eeglab redraw

% /////////////////////////////////////////////////////////////////////////
%%         ITC plots averaged over subjects for each channel
% /////////////////////////////////////////////////////////////////////////

CLim = [-.5 .5];

for i_event = 1:length(trialevent)
    
    for ch_itc = 1:length(electrode)

        % ITC values by electrode
        % (participants x conditions x events x electrodes x frequencies x timepoints)
        in_itc_chan = squeeze(mean(abs(itc(:,1,i_event,ch_itc,:,:)),1));
        out_itc_chan = squeeze(mean(abs(itc(:,2,i_event,ch_itc,:,:)),1));

        figure; 
        colormap('jet')

        % Subplot 1: out
        subplot(3,1,1); 
        imagesc(times,freqs,out_itc_chan,CLim);
        title({['ITC: ' char(trialevent{i_event}) ', Out']; char(elec_names(ch_itc))});  
        set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
        ylabel('Freq (Hz)'); 
        colorbar
        % Subplot 1: in
        subplot(3,1,2); 
        imagesc(times,freqs,in_itc_chan,CLim);
        title({['ITC: ' char(trialevent{i_event}) ', In']; char(elec_names(ch_itc))});
        set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
        ylabel('Freq (Hz)'); 
        colorbar
        % Subplot 1: out-in
        subplot(3,1,3); 
        imagesc(times,freqs,out_itc_chan-in_itc_chan,CLim); 
        title({['ITC: ' char(trialevent{i_event}) ', Out-In']; char(elec_names(ch_itc))});
        set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
        ylabel('Freq (Hz)'); xlabel('Time (ms)');
        colorbar

        %Overall subplot title
    %     supertitle(['ITC: ' char(elec_names(ch_itc))],'FontSize',12)

        clear in_itc_chan out_itc_chan 

    end
    clear ch_itc
end

% close all
% eeglab redraw
clear i_event

% .........................................................................   
% /////////////////////////////////////////////////////////////////////////
%%                            ERSP Top Plots
% /////////////////////////////////////////////////////////////////////////
% ......................................................................... 

%A Topoplot needs to collapse across frequency and time so it can show the 
% data across electrodes

% ......................................................................... 
% Set the range of frequency to consider
flims{1} = [1 3]; % delta
flims{2} = [4 7]; % theta
flims{3} = [8 12]; % alpha
flims{4} = [13 18]; % beta1
flims{5} = [19 22]; % beta2
flims{6} = [23 30]; % gamma

% flims{1} = [1 3]; % delta
% flims{1} = [4 7]; % theta
% flims{3} = [8 12]; % alpha
% flims{4} = [19 22]; % beta2
% flims{5} = [23 30]; % gamma

% ......................................................................... 
% Set the range of time to consider
% tlims{1} = [-400 -200];
% tlims{2} = [-200 0];
% tlims{3} = [0 200];
% tlims{4} = [200 400];
% tlims{5} = [400 600];

tlims{1} = [-200 -100];
tlims{2} = [-100 0];
tlims{3} = [0 100];
tlims{4} = [100 200];
tlims{5} = [200 300];
tlims{6} = [300 400];
tlims{7} = [400 500];
tlims{8} = [500 600];
tlims{9} = [600 700];

% ......................................................................... 
CLim = [-1.5 1.5];
% .........................................................................

for i_event = 1:length(trialevent) %loop through events

    for fq_i = 1:length(flims) %loop through set of frequencies

        iflims = flims{fq_i}; %select each frequency range
        %this code finds the frequencies you want from the freqs variable
        freq_lims = find(freqs>= iflims(1),1):find(freqs>= iflims(2),1)-1;
        
        if isempty(freq_lims) %in case the upper frequency < iflims(2)
            freq_lims = find(freqs>= iflims(1),1):find(freqs>= (iflims(2)-1),1)-1;
        end

        figure('OuterPosition',[313 537 1515 545]) %new figure for every frequency range

        for tl_i = 1:length(tlims) %loop through set of times

            itlims = tlims{tl_i}; %select each time range
            %this code finds the times you want from the timess variable
            time_lims = find(times>= itlims(1),1):find(times>= itlims(2),1)-1;

            % .....................................................................   
            %Here you need a 1D variable, electrodes.
            %By default it will take the mean across participants, events, times and frequencies, and show the data for each set
            % (participants x conditions x events x electrodes x frequencies x timepoints)
            in_ersp_top = squeeze(mean(mean(mean(ersp(:,1,i_event,electrode,freq_lims,...
                time_lims),5),6),1));
            out_ersp_top = squeeze(mean(mean(mean(ersp(:,2,i_event,electrode,freq_lims,...
                time_lims),5),6),1));
            diff_ersp_top = out_ersp_top - in_ersp_top;

            % .....................................................................
            % Creating plots

            % Subplot Out Condition
            subtightplot(3,length(tlims),tl_i,[0.01,0.01],[0.05,0.07],[0.05,0.05]);
%             subplot(3,length(tlims),tl_i);
            %This code creates the topoplots. You need to replace all the non-brain electrodes 
            % with NaN.
            topoplot([out_ersp_top' NaN NaN NaN],EEG.chanlocs,'maplimits',CLim,'plotrad',0.6,...
                'colormap',jet,'plotchans',electrode,'emarker',{'.','k',9,1}); colorbar;
            title(['Out: ' num2str(itlims(1)) ' to ' num2str(itlims(2)) ' ms'],'FontSize',8);

            % Subplot In Condition
            subtightplot(3,length(tlims),length(tlims)+tl_i,[0.01,0.01],[0.05,0.07],[0.05,0.05]);
            %This code creates the topoplots. You need to replace all the non-brain electrodes 
            % with NaN.
            topoplot([in_ersp_top' NaN NaN NaN],EEG.chanlocs,'maplimits',CLim,'plotrad',0.6,...
                'colormap',jet,'plotchans',electrode,'emarker',{'.','k',9,1}); colorbar;
            title(['In: ' num2str(itlims(1)) ' to ' num2str(itlims(2)) ' ms'],'FontSize',8);

            % Subplot Out-In Condition
            subtightplot(3,length(tlims),(length(tlims)*2)+tl_i,[0.01,0.01],[0.05,0.07],[0.05,0.05]);
            %This code creates the topoplots. You need to replace all the non-brain electrodes 
            % with NaN.
            topoplot([diff_ersp_top' NaN NaN NaN],EEG.chanlocs,'maplimits',CLim,'plotrad',0.6,...
                'colormap',jet,'plotchans',electrode,'emarker',{'.','k',9,1}); colorbar;
            title(['Out-In: ' num2str(itlims(1)) ' to ' num2str(itlims(2)) ' ms'],'FontSize',8);

            clear in_ersp_top out_ersp_top diff_ersp_top tflims time_lims itlims
        end

        % Overall subplot title
        supertitle({['ERSP: ' char(trialevent{i_event})]; [num2str(iflims(1)) ' to ' num2str(iflims(2)) ' Hz']},...
            'FontSize',12)

        clear iflims freq_lims
    end

end


clear tlims flims fq_i tl_i i_event

% close all
% eeglab redraw


% .........................................................................   
% /////////////////////////////////////////////////////////////////////////
%%                            ITC Top Plots
% /////////////////////////////////////////////////////////////////////////
% ......................................................................... 

%A Topoplot needs to collapse across frequency and time so it can show the 
% data across electrodes

% ......................................................................... 
% Set the range of frequency to consider
% flims{1} = [1 3]; % delta
% flims{1} = [4 7]; % theta
% flims{2} = [8 12]; % alpha
% flims{3} = [13 18]; % beta1
% flims{4} = [19 22]; % beta2
% flims{5} = [23 30]; % gamma

flims{1} = [1 3]; % delta
flims{2} = [4 7]; % theta
flims{3} = [8 12]; % alpha
flims{4} = [13 18]; % beta1

% ......................................................................... 
% Set the range of time to consider
tlims{1} = [0 100];
tlims{2} = [100 200];
tlims{3} = [200 300];
tlims{4} = [300 400];
tlims{5} = [400 500];

% tlims{1} = [0 50];
% tlims{2} = [50 100];
% tlims{3} = [100 150];
% tlims{4} = [150 200];
% tlims{5} = [200 250];
% tlims{6} = [250 300];
% tlims{7} = [300 350];
% tlims{8} = [350 400];

% ......................................................................... 
CLim = [-0.5 0.5];
% .........................................................................

for i_event = 1:length(trialevent) %loop through events
    
    for fq_i = 1:length(flims) %loop through set of frequencies

        iflims = flims{fq_i}; %select each frequency range
        %this code finds the frequencies you want from the freqs variable
        freq_lims = find(freqs>= iflims(1),1):find(freqs>= iflims(2),1)-1;
        
        if isempty(freq_lims) %in case the upper frequency < iflims(2)
            freq_lims = find(freqs>= iflims(1),1):find(freqs>= (iflims(2)-1),1)-1;
        end

        figure('OuterPosition',[313 537 1515 545]) %new figure for every frequency range

        for tl_i = 1:length(tlims) %loop through set of times

            itlims = tlims{tl_i}; %select each time range

            %this code finds the times you want from the timess variable
            time_lims = find(times>= itlims(1),1):find(times>= itlims(2),1)-1;

            % .....................................................................   
            %Here you need a 1D variable, electrodes.
            %By default it will take the mean across participants, events, times and frequencies, and show the data for each set
            % (participants x conditions x events x electrodes x frequencies x timepoints)
            in_itc_top = squeeze(mean(mean(mean(abs(itc(:,1,i_event,electrode,freq_lims,time_lims)),5),6),1));
            out_itc_top = squeeze(mean(mean(mean(abs(itc(:,2,i_event,electrode,freq_lims,time_lims)),5),6),1));
            diff_itc_top = out_itc_top - in_itc_top;

            % .....................................................................
            % Creating plots

            % Subplot Out Condition
            subtightplot(3,length(tlims),tl_i,[0.01,0.01],[0.05,0.07],[0.05,0.05]);
            %This code creates the topoplots. You need to replace all the non-brain electrodes 
            % with NaN.
            topoplot([out_itc_top' NaN NaN NaN],EEG.chanlocs,'maplimits',CLim,'plotrad',0.6,...
                'colormap',jet,'plotchans',electrode,'emarker',{'.','k',9,1}); colorbar;
            title(['Out: ' num2str(itlims(1)) ' to ' num2str(itlims(2)) ' ms'],'FontSize',8);

            % Subplot In Condition
            subtightplot(3,length(tlims),length(tlims)+tl_i,[0.01,0.01],[0.05,0.07],[0.05,0.05]);
            %This code creates the topoplots. You need to replace all the non-brain electrodes 
            % with NaN.
            topoplot([in_itc_top' NaN NaN NaN],EEG.chanlocs,'maplimits',CLim,'plotrad',0.6,...
                'colormap',jet,'plotchans',electrode,'emarker',{'.','k',9,1}); colorbar;
            title(['In: ' num2str(itlims(1)) ' to ' num2str(itlims(2)) ' ms'],'FontSize',8);

            % Subplot Out-In Condition
            subtightplot(3,length(tlims),(length(tlims)*2)+tl_i,[0.01,0.01],[0.05,0.07],[0.05,0.05]);
            %This code creates the topoplots. You need to replace all the non-brain electrodes 
            % with NaN.
            topoplot([diff_itc_top' NaN NaN NaN],EEG.chanlocs,'maplimits',CLim,'plotrad',0.6,...
                'colormap',jet,'plotchans',electrode,'emarker',{'.','k',9,1}); colorbar;
            title(['Out-In: ' num2str(itlims(1)) ' to ' num2str(itlims(2)) ' ms'],'FontSize',8);

            clear in_ersp_top out_ersp_top diff_ersp_top tflims time_lims itlims
        end

        % Overall subplot title
        supertitle({['ITC: ' char(trialevent{i_event})]; [num2str(iflims(1)) ' to ' num2str(iflims(2)) ' Hz']},...
            'FontSize',12)

        clear iflims freq_lims
    end

end

clear tlims flims fq_i tl_i i_event   
    
  

% /////////////////////////////////////////////////////////////////////////
%%  ERSP plots averaged over subjects and channels: Targets vs Standards
% /////////////////////////////////////////////////////////////////////////
   
% (participants x conditions x events x electrodes x frequencies x timepoints)
% standards
in_stand_ersp = squeeze(mean(mean(ersp(:,1,1,:,:,:),4),1));
out_stand_ersp = squeeze(mean(mean(ersp(:,2,1,:,:,:),4),1));
% targets
in_targ_ersp = squeeze(mean(mean(ersp(:,1,2,:,:,:),4),1));
out_targ_ersp = squeeze(mean(mean(ersp(:,2,2,:,:,:),4),1));
% targets - standards
in_edif_ersp = in_targ_ersp - in_stand_ersp;
out_edif_ersp = out_targ_ersp - out_stand_ersp;

figure; 
CLim = [-1.5 1.5];
colormap('jet')

% Subplot 1: out
subplot(3,1,1); 
imagesc(times,freqs,out_edif_ersp,CLim);
title('ERSP: Targets - Standards; Out'); 
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
ylabel('Freq (Hz)');
colorbar
% Subplot 1: in
subplot(3,1,2); 
imagesc(times,freqs,in_edif_ersp,CLim);
title('ERSP: Targets - Standards; In'); 
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
ylabel('Freq (Hz)');
colorbar
% Subplot 1: out-in
subplot(3,1,3); 
imagesc(times,freqs,out_edif_ersp-in_edif_ersp,CLim); 
title('ERSP: Targets - Standards; Out-In');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
ylabel('Freq (Hz)'); xlabel('Time (ms)');
colorbar

% Overall subplot title
% supertitle('ERSP: Targets - Standards','FontSize',12);
    
clear out_edif_ersp in_edif_ersp in_targ_ersp in_stand_ersp out_targ_ersp out_stand_ersp   

% /////////////////////////////////////////////////////////////////////////
%%   ITC plots averaged over subjects and channels: Targets vs Standards
% /////////////////////////////////////////////////////////////////////////


% (participants x conditions x events x electrodes x frequencies x timepoints)
% standards
in_stand_itc = squeeze(mean(mean(abs(itc(:,1,1,:,:,:)),4),1));
out_stand_itc = squeeze(mean(mean(abs(itc(:,2,1,:,:,:)),4),1));
% targets
in_targ_itc = squeeze(mean(mean(abs(itc(:,1,2,:,:,:)),4),1));
out_targ_itc = squeeze(mean(mean(abs(itc(:,2,2,:,:,:)),4),1));
% targets - standards
in_edif_itc = in_targ_itc - in_stand_itc;
out_edif_itc = out_targ_itc - out_stand_itc;


CLim = [-.5 .5];
figure; 
colormap('jet')

% Subplot 1: out
subplot(3,1,1); 
imagesc(times,freqs,out_edif_itc,CLim);
title('ITC: Targets - Standards; Out');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
ylabel('Freq (Hz)');
colorbar
% Subplot 1: in
subplot(3,1,2); 
imagesc(times,freqs,in_edif_itc,CLim);
title('ITC: Targets - Standards; In');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
ylabel('Freq (Hz)');
colorbar
% Subplot 1: out-in
subplot(3,1,3); 
imagesc(times,freqs,out_edif_itc-in_edif_itc,CLim); 
title('ITC: Targets - Standards; Out-In');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
ylabel('Freq (Hz)'); xlabel('Time (ms)');
colorbar

%Overall subplot title
% supertitle('ITC: Overall Average','FontSize',12)

clear out_edif_itc in_edif_itc in_targ_itc in_stand_itc out_targ_itc out_stand_itc 
    
% close all
% eeglab redraw    
    




%%


% #########################################################################
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% #########################################################################
%                            NORMALIZED PLOTS
% #########################################################################
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% #########################################################################




% .........................................................................   
% /////////////////////////////////////////////////////////////////////////
%%       Normalized ERSP plots averaged over subjects & channels
% /////////////////////////////////////////////////////////////////////////
% .........................................................................   


% ERSP values by condition and event
% (participants x conditions x events x electrodes x frequencies x timepoints)

% Standards
in_stand_ersp_temp = squeeze(mean(mean(ersp(:,1,1,:,:,:),4),1));
out_stand_ersp_temp = squeeze(mean(mean(ersp(:,2,1,:,:,:),4),1));
% Targets
in_targ_ersp_temp = squeeze(mean(mean(ersp(:,1,2,:,:,:),4),1));
out_targ_ersp_temp = squeeze(mean(mean(ersp(:,2,2,:,:,:),4),1));

%Here we are also going to take the difference from the average of the other conditions.
% diffs = [1:2]; diffs(1) = [];
diff_ersp_norm =  squeeze(mean(mean(mean(mean(ersp(:,:,:,:,:,:),1),2),3),4));

% .........................................................................  
% Create normalized standards and targets separately

% Normalize Standards
out_stand_ersp_norm = squeeze(out_stand_ersp_temp-diff_ersp_norm);
in_stand_ersp_norm = squeeze(in_stand_ersp_temp-diff_ersp_norm);

% Normalize Targets
out_targ_ersp_norm = squeeze(out_targ_ersp_temp-diff_ersp_norm);
in_targ_ersp_norm = squeeze(in_targ_ersp_temp-diff_ersp_norm);

% .........................................................................  
% Spectogram plot of normalized data
% .........................................................................  
%This variable sets the scale of the color axis, which corresponds to the itc or power values.
CLim = ([-1.5 1.5]);
% .........................................................................  

% Standards Figure
figure; 
colormap('jet')
% Subplot 1: out
subplot(3,1,1); 
imagesc(times,freqs,out_stand_ersp_norm,CLim);
title(['Normalized ERSP: Out; ' char(trialevent{1})]); 
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
ylabel('Freq (Hz)'); 
colorbar
% Subplot 2: in
subplot(3,1,2); 
imagesc(times,freqs,in_stand_ersp_norm,CLim);
title(['Normalized ERSP: In; ' char(trialevent{1})]); 
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
ylabel('Freq (Hz)'); 
colorbar
% Subplot 3: out-in
subplot(3,1,3);
imagesc(times,freqs,out_stand_ersp_norm-in_stand_ersp_norm,CLim); 
title(['Normalized ERSP: Out-In; ' char(trialevent{1})]);
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
ylabel('Freq (Hz)'); xlabel('Time (ms)');
colorbar


% Targets Figure
figure; 
colormap('jet')
% Subplot 1: out
subplot(3,1,1); 
imagesc(times,freqs,out_targ_ersp_norm,CLim);
title(['Normalized ERSP: Out; ' char(trialevent{2})]); 
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
ylabel('Freq (Hz)'); 
colorbar
% Subplot 2: in
subplot(3,1,2); 
imagesc(times,freqs,in_targ_ersp_norm,CLim);
title(['Normalized ERSP: In; ' char(trialevent{2})]); 
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
ylabel('Freq (Hz)'); 
colorbar
% Subplot 3: out-in
subplot(3,1,3);
imagesc(times,freqs,out_targ_ersp_norm-in_targ_ersp_norm,CLim); 
title(['Normalized ERSP: Out-In; ' char(trialevent{2})]);
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
ylabel('Freq (Hz)'); xlabel('Time (ms)');
colorbar


% Targets-Standards Figure
figure; 
colormap('jet')
% Subplot 1: out
subplot(3,1,1); 
imagesc(times,freqs,out_targ_ersp_norm-out_stand_ersp_norm,CLim);
title('Normalized ERSP: Out; Targets - Standards'); 
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
ylabel('Freq (Hz)'); 
colorbar
% Subplot 2: in
subplot(3,1,2); 
imagesc(times,freqs,in_targ_ersp_norm-in_stand_ersp_norm,CLim);
title('Normalized ERSP: In; Targets - Standards'); 
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
ylabel('Freq (Hz)'); 
colorbar
% Subplot 3: out-in
subplot(3,1,3);
imagesc(times,freqs,(out_targ_ersp_norm-out_stand_ersp_norm) - (in_targ_ersp_norm-in_stand_ersp_norm),...
    CLim); 
title('Normalized ERSP: Out-In; Targets - Standards');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
ylabel('Freq (Hz)'); xlabel('Time (ms)');
colorbar




clear diff_ersp_norm out_targ_ersp_norm in_targ_ersp_norm out_stand_ersp_norm...
    in_stand_ersp_norm in_stand_ersp_temp out_stand_ersp_temp in_targ_ersp_temp...
    out_targ_ersp_temp
  
% close all
% eeglab redraw




% .........................................................................   
% /////////////////////////////////////////////////////////////////////////
%%                        Normalized ERSP Top Plots
% /////////////////////////////////////////////////////////////////////////
% ......................................................................... 

%A Topoplot needs to collapse across frequency and time so it can show the 
% data across electrodes

% ......................................................................... 
% Set the range of frequency to consider
flims{1} = [1 3]; % delta
flims{2} = [4 7]; % theta
flims{3} = [8 12]; % alpha
flims{4} = [13 18]; % beta1
flims{5} = [19 22]; % beta2
flims{6} = [23 29]; % beta3
flims{7} = [30 35]; % gamma

% ......................................................................... 
% Set the range of time to consider
tlims{1} = [-300 -150];
tlims{2} = [-150 0];
tlims{3} = [0 150];
tlims{4} = [150 300];
tlims{5} = [300 450];
tlims{6} = [450 600];

% tlims{1} = [-100 0];
% tlims{2} = [0 100];
% tlims{3} = [100 200];
% tlims{4} = [200 300];
% tlims{5} = [300 400];
% tlims{6} = [400 500];
% tlims{7} = [500 600];
% tlims{8} = [600 700];

% ......................................................................... 
CLim = [-1.5 1.5];
% .........................................................................

% Create plots
for i_event = 1:length(trialevent) %loop through events

    for fq_i = 1:length(flims) %loop through set of frequencies

        iflims = flims{fq_i}; %select each frequency range
        %this code finds the frequencies you want from the freqs variable
        freq_lims = find(freqs>= iflims(1),1):find(freqs>= iflims(2),1)-1;
        
        if isempty(freq_lims) %in case the upper frequency < iflims(2)
            freq_lims = find(freqs>= iflims(1),1):find(freqs>= (iflims(2)-1),1)-1;
        end

        figure('OuterPosition',[311 397 1273 650]) %new figure for every frequency range

        for tl_i = 1:length(tlims) %loop through set of times

            itlims = tlims{tl_i}; %select each time range
            %this code finds the times you want from the timess variable
            time_lims = find(times>= itlims(1),1):find(times>= itlims(2),1)-1;
            
            % .....................................................................   
            % Difference ERSP for normalization
            diff_ersp_norm =  squeeze(mean(mean(mean(mean(mean(ersp(:,:,:,electrode,freq_lims,...
                time_lims),1),2),3),5),6));
            
            % .....................................................................   
            %By default it will take the mean across participants, events, times and frequencies, and show the data for each set
            % (participants x conditions x events x electrodes x frequencies x timepoints)
            in_normersp_top = squeeze(squeeze(mean(mean(mean(ersp(:,1,i_event,electrode,freq_lims,...
                time_lims),5),6),1)) - diff_ersp_norm);
            out_normersp_top = squeeze(squeeze(mean(mean(mean(ersp(:,2,i_event,electrode,freq_lims,...
                time_lims),5),6),1)) - diff_ersp_norm);
            diff_normersp_top = out_normersp_top - in_normersp_top;

            % .....................................................................
            % Creating plots

            % Subplot Out Condition
            subtightplot(3,length(tlims),tl_i,[0.01,0.01],[0.05,0.07],[0.05,0.05]);
%             subplot(3,length(tlims),tl_i);
            %This code creates the topoplots. You need to replace all the non-brain electrodes 
            % with NaN.
            topoplot([out_normersp_top' NaN NaN NaN],EEG.chanlocs,'maplimits',CLim,'plotrad',0.6,...
                'colormap',jet,'plotchans',electrode,'emarker',{'.','k',9,1}); colorbar;
            title(['Out: ' num2str(itlims(1)) ' to ' num2str(itlims(2)) ' ms'],'FontSize',8);

            % Subplot In Condition
            subtightplot(3,length(tlims),length(tlims)+tl_i,[0.01,0.01],[0.05,0.07],[0.05,0.05]);
            %This code creates the topoplots. You need to replace all the non-brain electrodes 
            % with NaN.
            topoplot([in_normersp_top' NaN NaN NaN],EEG.chanlocs,'maplimits',CLim,'plotrad',0.6,...
                'colormap',jet,'plotchans',electrode,'emarker',{'.','k',9,1}); colorbar;
            title(['In: ' num2str(itlims(1)) ' to ' num2str(itlims(2)) ' ms'],'FontSize',8);

            % Subplot Out-In Condition
            subtightplot(3,length(tlims),(length(tlims)*2)+tl_i,[0.01,0.01],[0.05,0.07],[0.05,0.05]);
            %This code creates the topoplots. You need to replace all the non-brain electrodes 
            % with NaN.
            topoplot([diff_normersp_top' NaN NaN NaN],EEG.chanlocs,'maplimits',CLim,'plotrad',0.6,...
                'colormap',jet,'plotchans',electrode,'emarker',{'.','k',9,1}); colorbar;
            title(['Out-In: ' num2str(itlims(1)) ' to ' num2str(itlims(2)) ' ms'],'FontSize',8);

            clear in_normersp_top out_normersp_top diff_normersp_top tflims time_lims itlims...
                diff_ersp_norm
        end

        % Overall subplot title
        supertitle({['Normalized ERSP: ' char(trialevent{i_event})]; [num2str(iflims(1)) ' to ' num2str(iflims(2)) ' Hz']},...
            'FontSize',12)
        

        clear iflims freq_lims
    end

end

clear tlims flims fq_i tl_i i_event

% close all
% eeglab redraw

% /////////////////////////////////////////////////////////////////////////
%%             Targets - Standards Normalized ERSP Top Plots
% /////////////////////////////////////////////////////////////////////////

%A Topoplot needs to collapse across frequency and time so it can show the 
% data across electrodes

% ......................................................................... 
% Set the range of frequency to consider
flims{1} = [1 3]; % delta
flims{2} = [4 7]; % theta
flims{3} = [8 12]; % alpha
flims{4} = [13 18]; % beta1
flims{5} = [19 22]; % beta2
flims{6} = [23 29]; % beta3
flims{7} = [30 35]; % gamma

% ......................................................................... 
% Set the range of time to consider
tlims{1} = [-300 -150];
tlims{2} = [-150 0];
tlims{3} = [0 150];
tlims{4} = [150 300];
tlims{5} = [300 450];
tlims{6} = [450 600];

% tlims{1} = [-100 0];
% tlims{2} = [0 100];
% tlims{3} = [100 200];
% tlims{4} = [200 300];
% tlims{5} = [300 400];
% tlims{6} = [400 500];
% tlims{7} = [500 600];
% tlims{8} = [600 700];

% ......................................................................... 
CLim = [-1.5 1.5];
% .........................................................................

% Plots
for fq_i = 1:length(flims) %loop through set of frequencies

    iflims = flims{fq_i}; %select each frequency range
    %this code finds the frequencies you want from the freqs variable
    freq_lims = find(freqs>= iflims(1),1):find(freqs>= iflims(2),1)-1;

    if isempty(freq_lims) %in case the upper frequency < iflims(2)
        freq_lims = find(freqs>= iflims(1),1):find(freqs>= (iflims(2)-1),1)-1;
    end

    figure('OuterPosition',[311 397 1273 650]) %new figure for every frequency range

    for tl_i = 1:length(tlims) %loop through set of times

        itlims = tlims{tl_i}; %select each time range
        %this code finds the times you want from the timess variable
        time_lims = find(times>= itlims(1),1):find(times>= itlims(2),1)-1;

        % .....................................................................   
        %By default it will take the mean across participants, events, times and frequencies, and show the data for each set
        % (participants x conditions x events x electrodes x frequencies x timepoints)
        % .....................................................................   
        % Difference ERSP for normalization
        diff_ersp_norm =  squeeze(mean(mean(mean(mean(mean(ersp(:,:,:,electrode,freq_lims,...
            time_lims),1),2),3),5),6));
        
        % .....................................................................
        % *Standards*
        in_stand_normersp_top = squeeze(squeeze(mean(mean(mean(ersp(:,1,1,electrode,freq_lims,...
            time_lims),5),6),1)) - diff_ersp_norm);
        out_stand_normersp_top = squeeze(squeeze(mean(mean(mean(ersp(:,2,1,electrode,freq_lims,...
            time_lims),5),6),1)) - diff_ersp_norm);
        diff_stand_normersp_top = out_stand_normersp_top - in_stand_normersp_top;
        
        % *Targets*
        in_targ_normersp_top = squeeze(squeeze(mean(mean(mean(ersp(:,1,2,electrode,freq_lims,...
            time_lims),5),6),1)) - diff_ersp_norm);
        out_targ_normersp_top = squeeze(squeeze(mean(mean(mean(ersp(:,2,2,electrode,freq_lims,...
            time_lims),5),6),1)) - diff_ersp_norm);
        diff_targ_normersp_top = out_targ_normersp_top - in_targ_normersp_top;

        % .....................................................................
        % Creating plots

        % Subplot Out Condition
        subtightplot(3,length(tlims),tl_i,[0.01,0.01],[0.05,0.07],[0.05,0.05]);
%             subplot(3,length(tlims),tl_i);
        %This code creates the topoplots. You need to replace all the non-brain electrodes 
        % with NaN.
        topoplot([(out_targ_normersp_top-out_stand_normersp_top)' NaN NaN NaN],...
            EEG.chanlocs,'maplimits',CLim,'plotrad',0.6,'colormap',jet,'plotchans',electrode,...
            'emarker',{'.','k',9,1}); colorbar;
        title(['Out: ' num2str(itlims(1)) ' to ' num2str(itlims(2)) ' ms'],'FontSize',9);

        % Subplot In Condition
        subtightplot(3,length(tlims),length(tlims)+tl_i,[0.01,0.01],[0.05,0.07],[0.05,0.05]);
        %This code creates the topoplots. You need to replace all the non-brain electrodes 
        % with NaN.
        topoplot([(in_targ_normersp_top-in_stand_normersp_top)' NaN NaN NaN],...
            EEG.chanlocs,'maplimits',CLim,'plotrad',0.6,'colormap',jet,'plotchans',electrode,...
            'emarker',{'.','k',9,1}); colorbar;
        title(['In: ' num2str(itlims(1)) ' to ' num2str(itlims(2)) ' ms'],'FontSize',9);

        % Subplot Out-In Condition
        subtightplot(3,length(tlims),(length(tlims)*2)+tl_i,[0.01,0.01],[0.05,0.07],[0.05,0.05]);
        %This code creates the topoplots. You need to replace all the non-brain electrodes 
        % with NaN.
        topoplot([(diff_targ_normersp_top - diff_stand_normersp_top)' NaN NaN NaN],EEG.chanlocs,...
            'maplimits',CLim,'plotrad',0.6,'colormap',jet,'plotchans',electrode,...
            'emarker',{'.','k',9,1}); colorbar;
        title(['Out-In: ' num2str(itlims(1)) ' to ' num2str(itlims(2)) ' ms'],'FontSize',9);

        clear tflims time_lims itlims in_targ_normersp_top in_stand_normersp_top...
            out_targ_normersp_top out_stand_normersp_top diff_stand_normersp_top...
            diff_targ_normersp_top diff_ersp_norm
    end

    % Overall subplot title
    supertitle({'Normalized ERSP: Targets - Standards'; [num2str(iflims(1)) ' to ' num2str(iflims(2)) ' Hz']},...
        'FontSize',12)


    clear iflims freq_lims
end


clear tlims flims fq_i tl_i i_event 


% close all
% eeglab redraw



% .........................................................................   
% /////////////////////////////////////////////////////////////////////////
%%       Normalized ERSP plots averaged over subjects for each channel
% /////////////////////////////////////////////////////////////////////////
% .........................................................................   


CLim = [-1.5 1.5];

for i_event = 1:length(trialevent) %loop through events

    for ch_ersp = 1:length(electrode)
        
        % .........................................................................  
        % ERSP values by electrode
        % (participants x conditions x events x electrodes x frequencies x timepoints)
        pref_ersp_chan = squeeze(mean(ersp(:,1,i_event,ch_ersp,:,:),1));
        npref_ersp_chan = squeeze(mean(ersp(:,2,i_event,ch_ersp,:,:),1));
        % .........................................................................  
        %Here we are also going to take the difference from the average of the other sets.
%         diff_ersp =  squeeze(mean(mean(mean(ersp(:,:,:,ch_ersp,:,:),1),2),3));
        diff_ersp =  squeeze(mean(mean(mean(mean(ersp(:,:,:,:,:,:),1),2),3),4));
        % .........................................................................  
        
        figure; 
        colormap('jet')

        % Subplot 1: out
        subplot(3,1,1); 
        imagesc(times,freqs,squeeze(npref_ersp_chan-diff_ersp),CLim);
        title({['Normalized ERSP: ' char(trialevent{i_event}) ', Out']; char(elec_names(ch_ersp))}); 
        set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
        ylabel('Freq (Hz)');
        colorbar
        % Subplot 1: in
        subplot(3,1,2); 
        imagesc(times,freqs,squeeze(pref_ersp_chan-diff_ersp),CLim);
        title({['Normalized ERSP: ' char(trialevent{i_event}) ', In']; char(elec_names(ch_ersp))}); 
        set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
        ylabel('Freq (Hz)');
        colorbar
        % Subplot 1: out-in
        subplot(3,1,3); 
        imagesc(times,freqs,(squeeze(npref_ersp_chan-diff_ersp)-squeeze(pref_ersp_chan-diff_ersp)),CLim); 
        title({['Normalized ERSP: ' char(trialevent{i_event}) ', Out - In']; char(elec_names(ch_ersp))});
        set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
        ylabel('Freq (Hz)'); xlabel('Time (ms)');
        colorbar
        
        clear in_ersp_chan out_ersp_chan diff_ersp

    end

    clear ch_ersp

end

clear i_event 

% close all
% eeglab redraw


% /////////////////////////////////////////////////////////////////////////
%                    Normalized: Targets - Standards
% /////////////////////////////////////////////////////////////////////////


CLim = [-1.5 1.5];

for ch_ersp = 1:length(electrode)
    
    % .........................................................................  
    % ERSP values by electrode
    % (participants x conditions x events x electrodes x frequencies x timepoints)
    % .........................................................................
    %Here we are also going to take the difference from the average of the other sets.
%     diff_ersp =  squeeze(mean(mean(mean(ersp(:,:,:,ch_ersp,:,:),2),3),1));
    diff_ersp =  squeeze(mean(mean(mean(mean(ersp(:,:,:,:,:,:),1),2),3),4));
    % .........................................................................  
    % Standards
    in_stand_ersp_chan = (squeeze(mean(mean(mean(ersp(:,1,1,ch_ersp,:,:),2),3),1))) - diff_ersp;
    out_stand_ersp_chan = (squeeze(mean(mean(mean(ersp(:,2,1,ch_ersp,:,:),2),3),1))) - diff_ersp;
    % .........................................................................  
    % Targets
    in_targ_ersp_chan = (squeeze(mean(mean(mean(ersp(:,1,2,ch_ersp,:,:),2),3),1))) - diff_ersp;
    out_targ_ersp_chan = (squeeze(mean(mean(mean(ersp(:,2,2,ch_ersp,:,:),2),3),1))) - diff_ersp;
    % .........................................................................
    % Targets - Standards
    in_diff_ersp_chan = in_targ_ersp_chan - in_stand_ersp_chan;
    out_diff_ersp_chan = out_targ_ersp_chan - out_stand_ersp_chan;
    % .........................................................................    
    
    figure; 
    colormap('jet')
    
    % Subplot 1: out
    subplot(3,1,1);
    imagesc(times,freqs,out_diff_ersp_chan,CLim);
    title({'Normalized ERSP: Targets - Standards';['Out, ' char(elec_names(ch_ersp))]}); 
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
    ylabel('Freq (Hz)');
    colorbar
    % Subplot 1: in
    subplot(3,1,2); 
    imagesc(times,freqs,in_diff_ersp_chan,CLim);
    title({'Normalized ERSP: Targets - Standards';['In, ' char(elec_names(ch_ersp))]});
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
    ylabel('Freq (Hz)');
    colorbar
    % Subplot 1: out-in
    subplot(3,1,3);  
    imagesc(times,freqs,out_diff_ersp_chan - in_diff_ersp_chan,CLim);
    title({'Normalized ERSP: Targets - Standards';['Out - In, ' char(elec_names(ch_ersp))]});
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
    ylabel('Freq (Hz)'); xlabel('Time (ms)');
    colorbar
    
    clear out_diff_ersp_chan in_diff_ersp_chan diff_ersp in_stand_ersp_chan...
        out_stand_ersp_chan in_targ_ersp_chan out_targ_ersp_chan

end
 
clear ch_ersp

% close all
% eeglab redraw


% .........................................................................   
% /////////////////////////////////////////////////////////////////////////
%%                       Normalized ITC Top Plots
% /////////////////////////////////////////////////////////////////////////
% ......................................................................... 

%A Topoplot needs to collapse across frequency and time so it can show the 
% data across electrodes

% ......................................................................... 
% Set the range of frequency to consider
% flims{1} = [1 3]; % delta
% flims{1} = [4 7]; % theta
% flims{2} = [8 12]; % alpha
% flims{3} = [13 18]; % beta1
% flims{4} = [19 22]; % beta2
% flims{5} = [23 30]; % gamma

flims{1} = [1 3]; % delta
flims{2} = [4 7]; % theta
flims{3} = [8 12]; % alpha
flims{4} = [13 18]; % beta1

% ......................................................................... 
% Set the range of time to consider
tlims{1} = [0 100];
tlims{2} = [100 200];
tlims{3} = [200 300];
tlims{4} = [300 400];
tlims{5} = [400 500];

% tlims{1} = [0 50];
% tlims{2} = [50 100];
% tlims{3} = [100 150];
% tlims{4} = [150 200];
% tlims{5} = [200 250];
% tlims{6} = [250 300];
% tlims{7} = [300 350];
% tlims{8} = [350 400];

% ......................................................................... 
CLim = [-0.5 0.5];
% .........................................................................

for i_event = 1:length(trialevent) %loop through events
    
    for fq_i = 1:length(flims) %loop through set of frequencies

        iflims = flims{fq_i}; %select each frequency range
        %this code finds the frequencies you want from the freqs variable
        freq_lims = find(freqs>= iflims(1),1):find(freqs>= iflims(2),1)-1;
        
        if isempty(freq_lims) %in case the upper frequency < iflims(2)
            freq_lims = find(freqs>= iflims(1),1):find(freqs>= (iflims(2)-1),1)-1;
        end

        figure('OuterPosition',[313 537 1515 545]) %new figure for every frequency range

        for tl_i = 1:length(tlims) %loop through set of times

            itlims = tlims{tl_i}; %select each time range

            %this code finds the times you want from the timess variable
            time_lims = find(times>= itlims(1),1):find(times>= itlims(2),1)-1;

            % .....................................................................   
            %Here you need a 1D variable, electrodes.
            %By default it will take the mean across participants, events, times and frequencies, and show the data for each set
            % (participants x conditions x events x electrodes x frequencies x timepoints)
            in_itc_top = squeeze(mean(mean(mean(abs(itc(:,1,i_event,electrode,freq_lims,time_lims)),5),6),1));
            out_itc_top = squeeze(mean(mean(mean(abs(itc(:,2,i_event,electrode,freq_lims,time_lims)),5),6),1));
            diff_itc_top = out_itc_top - in_itc_top;

            % .....................................................................
            % Creating plots

            % Subplot Out Condition
            subtightplot(3,length(tlims),tl_i,[0.01,0.01],[0.05,0.07],[0.05,0.05]);
            %This code creates the topoplots. You need to replace all the non-brain electrodes 
            % with NaN.
            topoplot([out_itc_top' NaN NaN NaN],EEG.chanlocs,'maplimits',CLim,'plotrad',0.6,...
                'colormap',jet,'plotchans',electrode,'emarker',{'.','k',9,1}); colorbar;
            title(['Out: ' num2str(itlims(1)) ' to ' num2str(itlims(2)) ' ms'],'FontSize',8);

            % Subplot In Condition
            subtightplot(3,length(tlims),length(tlims)+tl_i,[0.01,0.01],[0.05,0.07],[0.05,0.05]);
            %This code creates the topoplots. You need to replace all the non-brain electrodes 
            % with NaN.
            topoplot([in_itc_top' NaN NaN NaN],EEG.chanlocs,'maplimits',CLim,'plotrad',0.6,...
                'colormap',jet,'plotchans',electrode,'emarker',{'.','k',9,1}); colorbar;
            title(['In: ' num2str(itlims(1)) ' to ' num2str(itlims(2)) ' ms'],'FontSize',8);

            % Subplot Out-In Condition
            subtightplot(3,length(tlims),(length(tlims)*2)+tl_i,[0.01,0.01],[0.05,0.07],[0.05,0.05]);
            %This code creates the topoplots. You need to replace all the non-brain electrodes 
            % with NaN.
            topoplot([diff_itc_top' NaN NaN NaN],EEG.chanlocs,'maplimits',CLim,'plotrad',0.6,...
                'colormap',jet,'plotchans',electrode,'emarker',{'.','k',9,1}); colorbar;
            title(['Out-In: ' num2str(itlims(1)) ' to ' num2str(itlims(2)) ' ms'],'FontSize',8);

            clear in_ersp_top out_ersp_top diff_ersp_top tflims time_lims itlims
        end

        % Overall subplot title
        supertitle({['ITC: ' char(trialevent{i_event})]; [num2str(iflims(1)) ' to ' num2str(iflims(2)) ' Hz']},...
            'FontSize',12)

        clear iflims freq_lims
    end

end

clear tlims flims fq_i tl_i i_event   
    


    
    
    