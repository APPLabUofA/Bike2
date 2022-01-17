ccc
%
exp = 'Bike2';

subs = {'100' '101' '102'  '104' '106'  '108' '110'... 
        '114' '115' '116' '117' '118' '120' '121'...
        '122'  '126' '127' '129' '130' '131' '132' '133'...
         '135' '136'}; %new sample
     
nsubs = length(subs);
conds = {'sask' '110st' '83ave'};
conds_lab = {'Sask Drive Lane'; '110 Street Lane'; '83 Avenue Lane'};
nconds = length(conds);
Pathname = 'M:\Data\Bike_lanes\';
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
%%
for i_sub = 1:nsubs
    for i_cond = 1:nconds
        
        Filename = [subs{i_sub} '_' exp '_' conds{i_cond}];
        EEG = pop_loadset('filename',[Filename '_fft_Standard.set'],'filepath','M:\Data\Bike_lanes\segments_fft_JK\');
        %EEG = pop_loadset('filename',[Filename '_fft_Standard.set'],'filepath','M:\Data\Bike_lanes\segments_fft_JK_2\');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
    end
end
eeglab redraw
%%

 electrode = [1:16];

wavenumber = 6; %wavelet cycles 6 IS OG PARAMETER 
freqs = [1:.5:30]; %wavelet frequencies
power_out = [];
i_count = 0;
n_electrode = EEG.nbchan;
for i_sub = 1:nsubs
    fprintf(['Subject - ' num2str(i_sub) '. \n']);
    for i_cond = 1:nconds
        fprintf(['Condition - ' conds_lab{i_cond} '. \n']);
        i_count = i_count+1; %which data set in ALLEEG to use
        disp(ALLEEG(i_count).setname);
        power = [];
        n_trials = ALLEEG(i_count).trials;
        for i_trial = 1:n_trials
            for i_electrode = 1:length([electrode]) %pairs of contralateral electrodes
                tempdata = ALLEEG(i_count).data(electrode(i_electrode), :,i_trial) ;
                [temp_power temp_times temp_phase] = BOSC_tf(tempdata,freqs,EEG.srate,wavenumber);
                power(:,:,i_electrode,i_trial) = log(temp_power); %take log ?
            end
        end
        power_out(:,:,:,i_sub,i_cond) = mean(power,4); %save the power data, averaging across trials
    end
end

%% Collapse over time to get power spectra, average over subjects
power_spectra = squeeze(mean(power_out,2)); %collapse over time can pick time range here in 2nd dim.
%frequency x electrode x subject x condition

% power_spectra_diff = squeeze( power_spectra(:,left_electrode,:,:)-...
%                             power_spectra(:,right_electrode,:,:));
% power_spectra_diff_mean = squeeze(mean(power_spectra_diff,2));
% power_spectra_diff_se = squeeze(std(power_spectra_diff,[],2)/sqrt(nsubs));
power_spectra_mean = squeeze(mean(power_spectra,2));
power_spectra_se = squeeze(std(power_spectra,[],2)/sqrt(nsubs));


%% 
%Topographies

%%%now make topographies at each frequency bin
text_size = 8;
axis_label_size = 8;
title_size = 8;
line_width = 1;

F = 0:.5:30;
F = freqs;
delta_bins = find(F >= 0 & F <= 3);
theta_bins = find(F >= 4 & F <= 7);
alpha_bins = find(F >= 8 & F <= 12);
beta_bins = find(F >= 13 & F <= 30);

elec_locs = 'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced';

%OG CONDITIONS
%%%cond1%%%
topo_sask = squeeze(nanmean(nanmean(nanmean(power_out(alpha_bins,:,:,:,1),1),2),4));
topo_sask(17:18,1) = NaN;
%%%cond2%%%
topo_110st = squeeze(nanmean(nanmean(nanmean(power_out(alpha_bins,:,:,:,2),1),2),4));
topo_110st(17:18,1) = NaN;
%%%cond3%%%
topo_83ave = squeeze(nanmean(nanmean(nanmean(power_out(alpha_bins,:,:,:,3),1),2),4));
topo_83ave(17:18,1) = NaN;


plot_title = 'Power Topography';
%%%Difference%%%
%topo3 = topo1-topo2;
min_lim = min([topo_83ave;topo_83ave]);
max_lim = max([topo_83ave;topo_83ave]);
 diff_min_lim = min([topo_sask-topo_83ave]);
 diff_max_lim = max([topo_sask-topo_83ave]);


%original conds
figure;
subplot(1,3,1);
topoplot(topo_sask,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
    set(gca,'Color',[1 1 1]);
% title('Frequency Range 8-12 Hz');
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0]);
set(t,'LineWidth',line_width);

subplot(1,3,2);
topoplot(topo_110st,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
    set(gca,'Color',[1 1 1]);
%title(['Preferred CCW : ' plot_title]);
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0]);
set(t,'LineWidth',line_width);

subplot(1,3,3);
topoplot(topo_83ave,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
    set(gca,'Color',[1 1 1]);
%title(['Non-Preferred CW: ' plot_title]);
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0]);
%set(t,'LineWidth',line_width);


%IF MISSING ANYGHING, GO BACK TO SKATE TOPO JOINT ON M
