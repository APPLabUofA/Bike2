clear all
close all
%% ACCURACY
expacc = 'bike2';
subsacc = {'100' '101' '102'  '104' '106'  '108' '110'... 
        '114' '115' '116' '117' '118' '120' '121'...
        '122'  '126' '127' '129' '130' '131' '132' '133'...
         '135' '136'};
%%
nsubsacc = length(subsacc);
condsacc = {'sask'; '110st'; '83ave'};
new_condacc = {'Sask Drive'; '110 Street'; '83 Avenue'};

%preferred, clockwise - non-preffered, CCW
ncondsacc = length(condsacc);
Pathnameacc = 'M:\Data\Bike_lanes\';

% if ~exist([Pathname 'segments\'])
%     mkdir([Pathname 'segments\']);
% end
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% Marker Numbers
nStandardacc = 1;
nTargetacc = 2;
nFalseAlarmacc = 3;
nCorrectResponseacc = 4;

prop_correctacc = zeros(nsubsacc,length(new_condacc));
prop_correctRejacc = zeros(nsubsacc,length(new_condacc));
medianACC_correct = zeros(nsubsacc,length(new_condacc));
%medianRT_correct = zeros(nsubsacc,length(new_condacc));
%medianRT_falseAlarm = zeros(nsubsacc,length(new_condacc));
%%
for i_sub = 1:nsubsacc
    for i_cond = 1:ncondsacc
        
        Filenameacc = [subsacc{i_sub} '_' expacc '_' condsacc{i_cond} '.vhdr'];
        setnameacc = Filenameacc(1:end-5);
        
        EEG = pop_loadbv(Pathnameacc, Filenameacc);
        
        %% Find all event types and latencys
        
        %string trigger from brain recorder  (e.g. 'S  1')
        event_stringsacc = {EEG.event.type}; %array of marker label strings
        %time since recording start in ms (integer)
        event_latencyacc = [EEG.event.latency];
                
        %convert strings to integers so they are easier to work with
        event_markersacc = zeros(size(event_stringsacc));
        event_markersacc(strcmp(event_stringsacc,'S  1')) = nStandardacc;   %standard
        event_markersacc(strcmp(event_stringsacc,'S  2')) = nTargetacc;   %target
        event_markersacc(strcmp(event_stringsacc,'S  3')) = nFalseAlarmacc;   %false alarm
        event_markersacc(strcmp(event_stringsacc,'S  4')) = nCorrectResponseacc;   %correct response
        
       
        event_latencyacc(event_markersacc == 0) = []; %remove any extra triggers
        event_markersacc(event_markersacc == 0) = []; % remove any extra triggers
        
        %% now step through the arrays and check for stuff
        
        %setup counters
        count_tonesacc = 0;
        count_targetsacc = 0;
        count_standardsacc = 0;
        
        count_correctacc = 0;
        count_missesacc = 0;
        count_correctRejacc = 0;
        count_falseAlarmacc = 0;
        
        ACC_correct = [];
%         RT_correct = [];
%         RT_falseAlarm = [];
        
        %for every event
        for i_event = 1:length(event_markersacc)-1 %last one is a filler markers
            
            this_markeracc = event_markersacc(i_event);
            tone_timeacc = event_latencyacc(i_event);
            next_markeracc = event_markersacc(i_event+1);
            next_timeacc = event_latencyacc(i_event+1);
            potential_ACC = count_correctacc/50*100;
%             potential_RT = next_time-tone_time;
            
            %if it is a tone (|| means or)
            if this_markeracc == nTargetacc || this_markeracc == nStandardacc
                count_tonesacc = count_tonesacc + 1;
                fprintf('\n Tone Number: ') %\n is a new line
                fprintf(num2str(count_tonesacc))
                fprintf(' --> ')
                
                fprintf('This marker: ')
                fprintf(num2str(this_markeracc))
                fprintf(' , ')
                
                fprintf('Next marker: ')
                fprintf(num2str(next_markeracc))
                fprintf(' , ')
                
                
            end
            
            if this_markeracc == nTargetacc 
                count_targetsacc = count_targetsacc + 1;
                
                
                %if correct response
                if next_markeracc == nCorrectResponseacc 
                    count_correctacc = count_correctacc + 1;
                    ACC_correct = [ACC_correct potential_ACC];
                    fprintf('Responded -- > Cummulative Accuracy = ')
                    fprintf(num2str(potential_ACC))
                    fprintf(' %')
                    
                    %if miss since next is another tone
                elseif next_markeracc == nStandardacc || next_markeracc == nTargetacc
                    count_missesacc = count_missesacc + 1;
                    fprintf('Did not respond')
                    
                    %anything else?
                else
                    fprintf('Not 9 or 3 or 5 (or 4, 1 or 2)')
                    
                end
                
            elseif this_markeracc == nStandardacc
                count_standardsacc = count_standardsacc + 1;
                
                %if correct rejection since next is another tone
                if next_markeracc == nStandardacc || next_markeracc == nTargetacc
                    count_correctRejacc = count_correctRejacc + 1;
                    fprintf('Correct Rejection')
                    
                    %if false alarm
                elseif next_markeracc == nFalseAlarmacc
                    %RT_falseAlarm = [RT_falseAlarm potential_RT];
                    count_falseAlarmacc = count_falseAlarmacc + 1;
                    fprintf('False Alarm')

                    
                    %anything else?
                else
                    fprintf('Not 7 or 3 or 5 (or 3 1 2')
                end
                
            end %if target,elseStandard
        end %every event
        
                prop_correctacc(i_sub,i_cond) = count_correctacc / count_targetsacc;
        prop_correctRejacc(i_sub,i_cond) = count_correctRejacc / count_standardsacc;
        medianACC_correct(i_sub,i_cond) = median(potential_ACC);
        %medianRT_falseAlarm(i_sub,new_cond_index) = median(RT_falseAlarm);
        
    end
end

%%
n_conditionsacc = size(medianACC_correct,2);
%this the grand mean over subjects of the median RTs
grand_mean_ACC_Corr = mean(medianACC_correct);
%these are normal error bars (Standard Error)
grand_SE_ACC_Corr = std(medianACC_correct)/sqrt(nsubsacc);
%these are smaller within subject error bars
%made by subtracting away each subjects average
%from their other scores to remove between subject difference
sub_mean_ACC_Corr = mean(medianACC_correct,2); %average for each subject
%this subtracts each subjects average from their scores
%repmat repeats the matrix 4 times for each condition
mean_ACC_Corr_deviation = medianACC_correct - repmat(sub_mean_ACC_Corr,1,n_conditionsacc);

%calculating standard errors for plots
SE_Pacc = nanstd(mean(medianACC_correct(:,1:2),2))/sqrt(nsubsacc); 
SE_NPacc = nanstd(mean(medianACC_correct(:,3:4),2))/sqrt(nsubsacc);

%then take the standard error of those deviatoins from the mean
grand_withinSE_ACC_Corr = std(mean_ACC_Corr_deviation)/sqrt(nsubsacc);

%now do the same for proportion correct
grand_mean_prop_corracc = nanmean(prop_correctacc);
grand_SE_prop_corracc = nanstd(prop_correctacc)/sqrt(nsubsacc);
sub_mean_prop_corracc = nanmean(prop_correctacc,2);
prop_corr_deviationacc = prop_correctacc - repmat(sub_mean_prop_corracc,1,n_conditionsacc);
grand_withinSE_prop_corr = nanstd(prop_corr_deviationacc)/sqrt(nsubsacc);

%% original 4 conditions plus proportion correct
% %plot it
conds_plot = {'Saskatchewan Lane'; '110 St. Lane';'83 Ave. Lane'}; 
figure;
subplot (1,2,1)
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(grand_mean_ACC_Corr,grand_withinSE_ACC_Corr);
 ylim([90 100])
ylabel('Mean Accuracy')
xlabel('Condition')
title('Target Accuracy')
legend(conds_plot)
subplot(1,2,2)
barweb(grand_mean_prop_corracc,grand_withinSE_prop_corr);
 ylim([.9 1])
ylabel('Proportion')
xlabel('Condition')
title('Proportion of Targets responded to')

here = [100, 200,300];
there = [10, 20, 30];
barweb(100,10);


%Jasp variables for excel
acc_sask = medianACC_correct(:,1);
acc_110 = medianACC_correct(:,2);
acc_83 = medianACC_correct(:,3);



%% REACTION TIME

%%
ccc
exp = 'bike2';
subs = {'100' '101' '102'  '104' '106'  '108' '110'... 
        '114' '115' '116' '117' '118' '120' '121'...
        '122'  '126' '127' '129' '130' '131' '132' '133'...
         '135' '136'};

%%
nsubs = length(subs);
conds = {'sask'; '110st'; '83ave'};
new_conds = {'Sask Drive'; '110 Street'; '83 Avenue'};

%preferred, clockwise - non-preffered, CCW
nconds = length(conds);
Pathname = 'M:\Data\Bike_Mazumder\';

% if ~exist([Pathname 'segments\'])
%     mkdir([Pathname 'segments\']);
% end
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% Marker Numbers
nStandard = 1;
nTarget = 2;
nFalseAlarm = 3;
nCorrectResponse = 4;

prop_correct = zeros(nsubs,length(new_conds));
prop_correctRej = zeros(nsubs,length(new_conds));
medianRT_correct = zeros(nsubs,length(new_conds));
medianRT_falseAlarm = zeros(nsubs,length(new_conds));

%%
for i_sub = 1:nsubs
    for i_cond = 1:nconds
        
        Filename = [subs{i_sub} '_' exp '_' conds{i_cond} '.vhdr'];
        setname = Filename(1:end-5);
        
        EEG = pop_loadbv(Pathname, Filename);
        
        %% Find all event types and latencys
        
        %string trigger from brain recorder  (e.g. 'S  1')
        event_strings = {EEG.event.type}; %array of marker label strings
        %time since recording start in ms (integer)
        event_latency = [EEG.event.latency];
                
        %convert strings to integers so they are easier to work with
        event_markers = zeros(size(event_strings));
        event_markers(strcmp(event_strings,'S  1')) = nStandard;   %standard
        event_markers(strcmp(event_strings,'S  2')) = nTarget;   %target
        event_markers(strcmp(event_strings,'S  3')) = nFalseAlarm;   %false alarm
        event_markers(strcmp(event_strings,'S  4')) = nCorrectResponse;   %correct response
                  
        event_latency(event_markers == 0) = []; %remove any extra triggers
        event_markers(event_markers == 0) = []; % remove any extra triggers
        
        %% now step through the arrays and check for stuff
        
        %setup counters
        count_tones = 0;
        count_targets = 0;
        count_standards = 0;
        
        count_correct = 0;
        count_misses = 0;
        count_correctRej = 0;
        count_falseAlarm = 0;
        
        RT_correct = [];
        RT_falseAlarm = [];
        
        %for every event
        for i_event = 1:length(event_markers)-1 %last one is a filler markers
            
            this_marker = event_markers(i_event);
            tone_time = event_latency(i_event);
            next_marker = event_markers(i_event+1);
            next_time = event_latency(i_event+1);
            potential_RT = next_time-tone_time;
            
            %if it is a tone (|| means or)
            if this_marker == nTarget || this_marker == nStandard
                count_tones = count_tones + 1;
                fprintf('\n Tone Number: ') %\n is a new line
                fprintf(num2str(count_tones))
                fprintf(' --> ')
                
                fprintf('This marker: ')
                fprintf(num2str(this_marker))
                fprintf(' , ')
                
                fprintf('Next marker: ')
                fprintf(num2str(next_marker))
                fprintf(' , ')
                
                
            end
            
            if this_marker == nTarget
                count_targets = count_targets + 1;
                
                
                %if correct response
                if next_marker == nCorrectResponse 
                    count_correct = count_correct + 1;
                    RT_correct = [RT_correct potential_RT];
                    fprintf('Responded -- > RT = ')
                    fprintf(num2str(potential_RT))
                    fprintf(' ms')
                    
                    %if miss since next is another tone
                elseif next_marker == nStandard || next_marker == nTarget
                    count_misses = count_misses + 1;
                    fprintf('Did not respond')
                    
                    %anything else?
                else
                    fprintf('Not 9 or 3 or 5 (or 4, 1 or 2)')
                    
                end
                
            elseif this_marker == nStandard
                count_standards = count_standards + 1;
                
                %if correct rejection since next is another tone
                if next_marker == nStandard || next_marker == nTarget
                    count_correctRej = count_correctRej + 1;
                    fprintf('Correct Rejection')
                    
                    %if false alarm
                elseif next_marker == nFalseAlarm
                    RT_falseAlarm = [RT_falseAlarm potential_RT];
                    count_falseAlarm = count_falseAlarm + 1;
                    fprintf('False Alarm -- > RT = ')
                    fprintf(num2str(potential_RT))
                    fprintf(' ms')
                    
                    %anything else?
                else
                    fprintf('Not 7 or 3 or 5 (or 3 1 2')
                end
                
            end %if target,elseStandard
        end %every event
        
        prop_correct(i_sub,i_cond) = count_correct / count_targets;
        prop_correctRej(i_sub,i_cond) = count_correctRej / count_standards;
        medianRT_correct(i_sub,i_cond) = median(RT_correct);
        medianRT_falseAlarm(i_sub,i_cond) = median(RT_falseAlarm);
        
    end
end

%%

n_conditions = size(medianRT_correct,2);
%this the grand mean over subjects of the median RTs
grand_mean_RT_Corr = mean(medianRT_correct);
%these are normal error bars (Standard Error)
grand_SE_RT_Corr = std(medianRT_correct)/sqrt(nsubs);
%these are smaller within subject error bars
%made by subtracting away each subjects average
%from their other scores to remove between subject difference
sub_mean_RT_Corr = mean(medianRT_correct,2); %average for each subject
%this subtracts each subjects average from their scores
%repmat repeats the matrix 4 times for each condition
mean_RT_Corr_deviation = medianRT_correct - repmat(sub_mean_RT_Corr,1,n_conditions);
%then take the standard error of those deviatoins from the mean
grand_withinSE_RT_Corr = std(mean_RT_Corr_deviation)/sqrt(nsubs);


%now do the same for proportion correct
grand_mean_prop_corr = mean(prop_correct);
grand_SE_prop_corr = std(prop_correct)/sqrt(nsubs);
sub_mean_prop_corr = mean(prop_correct,2);
prop_corr_deviation = prop_correct - repmat(sub_mean_prop_corr,1,n_conditions);
grand_withinSE_prop_corr = std(prop_corr_deviation)/sqrt(nsubs);
%% original 4 conditions plus proportion correct
%plot it
%conds_plot = {'Sask'; '110 Street';'83 Ave'}; 
figure;
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
subplot(1,2,1)
barweb(grand_mean_RT_Corr,grand_withinSE_RT_Corr);
ylim([450 650])
ylabel('Median RT (ms)')
xlabel('Condition')
title('Target Reaction Time (w/i subject SE)')
legend(conds_plot)
subplot(1,2,2)
barweb(grand_mean_prop_corr,grand_withinSE_prop_corr);
ylim([.9 1])
ylabel('Proportion')
xlabel('Condition')
title('Proportion of Targets responded to')

%Jasp variables for excel
rt_sask = medianRT_correct(:,1);
rt_110 = medianRT_correct(:,2);
rt_83 = medianRT_correct(:,3);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 3 MANUSCRIPT 
%% Scatter plots for ACC

%%% Set some variables for line width and text size.
text_size = 12;
axis_label_size = 14;
title_size = 16;
legend_size = 10;
line_width = 3;

%%% Create our figure.
figure; hold on;

%%% Specify axis limits and labels.
xlim([0,5]);
ylim([40,100]);
xticks([0, 1, 2, 3, 4, 5]);
xticklabels([]);
title('Median Target Accuracy','FontSize', title_size,'FontWeight', 'bold');
xlabel('Condition','FontSize', axis_label_size,'FontWeight', 'bold');
ylabel('Accuracy','FontSize', axis_label_size,'FontWeight', 'bold');


%%% Change line width and other plot properties.
set(gca, 'linewidth', line_width, 'FontSize', text_size,...
    'FontWeight', 'bold', 'linewidth', line_width, 'box', 'off',...
    'color', 'none', 'Layer', 'Top');


%%% Determine colours for each condition.
colours = {'b', 'r', 'g', 'm'};

%%% Now loop through our conditions.
for i_plot = 1:4
    
    %%% Determine how many participants we have for this condition.
    n_parts = ones(length(medianACC_correct(:,1)), 1) * i_plot;
    
    %%% Create our scatter plot.
    h(i_plot) = scatter(n_parts, medianACC_correct(:,i_plot), 75, 'o', ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', colours{i_plot},...
        'LineWidth', 1.5);
    
    %%% Now plot the median of our particiapnt data.
    scatter(i_plot, median(medianACC_correct(:,i_plot)), 75, 'd', ...
        'MarkerEdgeColor', [0.5, 0.5, 0.5], 'MarkerFaceColor', 'y',...
        'LineWidth', 1.5);

end

%%% Make a legend for the participant data only.
legend(h, {'Preferred Clockwise', 'Preferred Counter-Clockwise', ...
    'Non-preferred Clockwise', 'Non-preferred Counter-Clockwise'},...
    'FontSize', legend_size);
legend('boxoff');

hold off;


%% Scatter plots for RT

%%% Set some variables for line width and text size.
text_size = 12;
axis_label_size = 14;
title_size = 16;
legend_size = 10
line_width = 3;

%%% Create our figure.
figure; hold on;

%%% Specify axis limits and labels.
xlim([0,5]);
% ylim([40,100]);
xticks([0, 1, 2, 3, 4, 5]);
xticklabels([]);
title('Median Target Response Time','FontSize', title_size,'FontWeight', 'bold');
xlabel('Condition','FontSize', axis_label_size,'FontWeight', 'bold');
ylabel('Response Time (ms)','FontSize', axis_label_size,'FontWeight', 'bold');


%%% Change line width and other plot properties.
set(gca, 'linewidth', line_width, 'FontSize', text_size,...
    'FontWeight', 'bold', 'linewidth', line_width, 'box', 'off',...
    'color', 'none', 'Layer', 'Top');


%%% Determine colours for each condition.
colours = {'b', 'r', 'g', 'm'};

%%% Now loop through our conditions.
for i_plot = 1:4
    
    %%% Determine how many participants we have for this condition.
    n_parts = ones(length(medianRT_correct(:,1)), 1) * i_plot;
    
    %%% Create our scatter plot.
    h(i_plot) = scatter(n_parts, medianRT_correct(:,i_plot), 75, 'o', ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', colours{i_plot},...
        'LineWidth', 1.5);
    
    %%% Now plot the median of our particiapnt data.
    scatter(i_plot, median(medianRT_correct(:,i_plot)), 75, 'd', ...
        'MarkerEdgeColor', [0.5, 0.5, 0.5], 'MarkerFaceColor', 'y',...
        'LineWidth', 1.5);

end

%%% Make a legend for the participant data only.
legend(h, {'Preferred Clockwise', 'Preferred Counter-Clockwise', ...
    'Non-preferred Clockwise', 'Non-preferred Counter-Clockwise'},...
    'FontSize', legend_size);
legend('boxoff');

hold off;

%%
%%%%%%%%%%%%%%%%%%%%%%
%%%% by averaged preference

%NEW PLOTS revision
%median accuracy by preference 

pref_acc = mean(medianACC_correct(:,1:2),2);
non_pref_acc = mean(medianACC_correct(:,3:4),2);
medianACC_globalpref = [pref_acc non_pref_acc];


%%% Set some variables for line width and text size.
text_size = 12;
axis_label_size = 14;
title_size = 16;
legend_size = 10;
line_width = 3;

%%% Create our figure.
figure; hold on;

%%% Specify axis limits and labels.
xlim([0, 3]);
ylim([40,100]);
xticks([1, 2]);
xticklabels([]);
title('Accuracy By Preference','FontSize', title_size,'FontWeight', 'bold');
xlabel('Condition','FontSize', axis_label_size,'FontWeight', 'bold');
ylabel('Accuracy','FontSize', axis_label_size,'FontWeight', 'bold');


%%% Change line width and other plot properties.
set(gca, 'linewidth', line_width, 'FontSize', text_size,...
    'FontWeight', 'bold', 'linewidth', line_width, 'box', 'off',...
    'color', 'none', 'Layer', 'Top');


%%% Determine colours for each condition.
colours = {'b', 'r'};

%%% Now loop through our conditions.
for i_plot = 1:2
    
    %%% Determine how many participants we have for this condition.
    n_parts = ones(length(medianACC_globalpref(:,1)), 1) * i_plot;
    
    %%% Create our scatter plot.
    h(i_plot) = scatter(n_parts, medianACC_globalpref(:,i_plot), 75, 'o', ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', colours{i_plot},...
        'LineWidth', 1.5);
    
    %%% Now plot the median of our particiapnt data.
    scatter(i_plot, median(medianACC_globalpref(:,i_plot)), 75, 'd', ...
        'MarkerEdgeColor', [0.5, 0.5, 0.5], 'MarkerFaceColor', 'y',...
        'LineWidth', 1.5);

end

legend(h, {'Preferred', 'Non-preferred'}, 'FontSize', legend_size);
legend('boxoff');

hold off;

%% Reaction time by preference
pref_RT = mean(medianRT_correct(:,1:2),2);
non_pref_rt = mean(medianRT_correct(:,3:4),2);
medianRT_global = [pref_RT non_pref_rt];

%%% Set some variables for line width and text size.
text_size = 12;
axis_label_size = 14;
title_size = 16;
legend_size = 10;
line_width = 3;

%%% Create our figure.
figure; hold on;

%%% Specify axis limits and labels.
xlim([0,3]);
% ylim([40,100]);
xticks([1, 2]);
xticklabels([]);
title('Response Time by Preference','FontSize', title_size,'FontWeight', 'bold');
xlabel('Condition','FontSize', axis_label_size,'FontWeight', 'bold');
ylabel('Response Time (ms)','FontSize', axis_label_size,'FontWeight', 'bold');


%%% Change line width and other plot properties.
set(gca, 'linewidth', line_width, 'FontSize', text_size,...
    'FontWeight', 'bold', 'linewidth', line_width, 'box', 'off',...
    'color', 'none', 'Layer', 'Top');


%%% Determine colours for each condition.
colours = {'b', 'r'};

%%% Now loop through our conditions.
for i_plot = 1:2
    
    %%% Determine how many participants we have for this condition.
    n_parts = ones(length(medianRT_global(:,1)), 1) * i_plot;
    
    %%% Create our scatter plot.
    h(i_plot) = scatter(n_parts, medianRT_global(:,i_plot), 75, 'o', ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', colours{i_plot},...
        'LineWidth', 1.5);
    
    %%% Now plot the median of our particiapnt data.
    scatter(i_plot, median(medianRT_correct(:,i_plot)), 75, 'd', ...
        'MarkerEdgeColor', [0.5, 0.5, 0.5], 'MarkerFaceColor', 'y',...
        'LineWidth', 1.5);

end

%%% Make a legend for the participant data only.
legend(h, {'Preferred', 'Non-preferred'},...
    'FontSize', legend_size);
legend('boxoff');

hold off;

