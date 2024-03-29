clear all
close all
ccc

exp = 'Bike2';
subs = {'100' '101' '102'  '104' '106'  '108' '110'... 
        '114' '115' '116' '117' '118' '120' '121'...
        '122'  '126' '127' '129' '130' '131' '132' '133'...
         '135' '136'}; %new sample
     
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
%         EEG = pop_loadset('filename',[Filename '_Targets.set'],'filepath','M:\Data\Bike_lanes\segments_JK_2\');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        EEG = pop_loadset('filename',[Filename '_Standards.set'],'filepath','M:\Data\Bike_lanes\segments_JK');
%         EEG = pop_loadset('filename',[Filename '_Standards.set'],'filepath','M:\Data\Bike_lanes\segments_JK_2');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
    end
end
eeglab redraw
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%subject erps
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

erp_out = [];
for i_sub = 1:nsubs
    fprintf(['Subject ' num2str(i_sub)])
    for i_cond = 1:nconds
        %average over trials (3rd dimension)
        erp_out(:,1,:,i_cond,i_sub) = mean(ALLEEG(1+ 2*((i_sub-1)*nconds+(i_cond-1))).data,3)'; %Targets
        erp_out(:,2,:,i_cond,i_sub) = mean(ALLEEG(2+ 2*((i_sub-1)*nconds+(i_cond-1))).data,3)'; %standards
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%grand average plots + difference waves
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

electrode = 15;
erp_diff_out = squeeze(erp_out(:,1,:,:,:)-erp_out(:,2,:,:,:));
figure('Color',[1 1 1]);
for i_cond = 1:nconds
    switch i_cond
        case 1
            colour = 'b';
        case 2
            colour = 'g';
        case 3
            colour = 'r';
            %         case 4
            %             colour = 'm';
    end
    
    subplot(2,nconds,i_cond);
    boundedline(EEG.times,squeeze(mean(erp_out(:,1,electrode,i_cond,:),5)),squeeze(std(erp_out(:,1,electrode,i_cond,:),[],5))./sqrt(nsubs),colour,...
        EEG.times,squeeze(mean(erp_out(:,2,electrode,i_cond,:),5)),squeeze(std(erp_out(:,2,electrode,i_cond,:),[],5))./sqrt(nsubs),'k');
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    if i_cond == 2
        legend('Targets','Standards','Location','NorthEast','Autoupdate', 'off');
    end
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title(conds{i_cond});
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    
    subplot(2,nconds,nconds+i_cond);
    boundedline(EEG.times,squeeze(mean(erp_diff_out(:,electrode,i_cond,:),4)),squeeze(std(erp_diff_out(:,electrode,i_cond,:),[],4))./sqrt(nsubs),colour);
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    if i_cond == 2
        legend('Targets-Standards','Location','NorthEast','Autoupdate', 'off');
    end
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title(conds{i_cond});
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% difference waves on same axis
%
% MANUSCRIPT: DIFFERENCE-WAVES P3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

electrode = 15;
erp_diff_out = squeeze(erp_out(:,1,:,:,:)-erp_out(:,2,:,:,:));

figure
boundedline(EEG.times,squeeze(mean(erp_diff_out(:,electrode,1,:),4)),squeeze(std(erp_diff_out(:,electrode,1,:),[],4))./sqrt(nsubs),'b',...%erp_diff_out(:,electrode,1,:),[],4
    EEG.times,squeeze(mean(erp_diff_out(:,electrode,2,:),4)),squeeze(std(erp_diff_out(:,electrode,2,:),[],4))./sqrt(nsubs),'g',...
    EEG.times,squeeze(mean(erp_diff_out(:,electrode,3,:),4)),squeeze(std(erp_diff_out(:,electrode,3,:),[],4))./sqrt(nsubs),'r');

set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Difference-Wave');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
L(1) = plot(nan, nan, 'b');
L(2) = plot(nan, nan, 'g');
L(3) = plot(nan, nan, 'r');
legend(L, {'Sask Drive', '110 Street', '83 Avenue'},'Location','northwest', 'Autoupdate', 'off')
% fill([300;300;450;450],[-8;12;12;-8],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
fill([355;355;505;505],[-8;12;12;-8],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
clear electrode;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                  P3 DIFFERENCE WAVE STATS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TARGET VARIABLES FOR JASP
electrode = 15;
erp_diff_out = squeeze(erp_out(:,1,:,:,:)-erp_out(:,2,:,:,:));
time_window = find(EEG.times>355,1)-1:find(EEG.times>505,1)-2;
P3_diffs = [squeeze(mean(erp_diff_out(time_window,electrode,1,:),1)),...
    squeeze(mean(erp_diff_out(time_window,electrode,2,:),1)),...
    squeeze(mean(erp_diff_out(time_window,electrode,3,:),1))];
writematrix(P3_diffs, 'Pz_P300_dif-wave.csv')
clear time_window;

%%
%t-test for dif waves
time_window = find(EEG.times>325,1)-1:find(EEG.times>450,1)-2;
erp_diff_out = squeeze(erp_out(:,1,:,:,:)-erp_out(:,2,:,:,:));

%%%timepoints X events X electrodes X conditions X participants%%%%

%  COMPARING DIFFERENCE WAVES for non-preferred conditions
%[h p ci stat] = ttest(squeeze(mean(erp_diff_out(time_window,electrode,1,:),1)),squeeze(mean(erp_diff_out(time_window,electrode,2,:),1)),.05,'both',1)

%non-parametric
%[p, h,stats] = ranksum(squeeze(mean(erp_diff_out(time_window,electrode,2,:),1)),squeeze(mean(erp_diff_out(time_window,electrode,1,:),1)),'alpha', 0.05, 'tail','both')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         TARGET ERP FIGURES + STATS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Targets on same axis
electrode = 13;

figure;
boundedline(EEG.times,squeeze(mean(erp_out(:,1,electrode,1,:),5)), squeeze(std(erp_out(:,1,electrode,1,:),[],5))./sqrt(nsubs),'b',...
    EEG.times,squeeze(mean(erp_out(:,1,electrode,2,:),5)), squeeze(std(erp_out(:,1,electrode,2,:),[],5))./sqrt(nsubs),'g',...
    EEG.times,squeeze(mean(erp_out(:,1,electrode,3,:),5)), squeeze(std(erp_out(:,1,electrode,3,:),[],5))./sqrt(nsubs),'r')

set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Target Tones');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
L(1) = plot(nan, nan, 'b');
L(2) = plot(nan, nan, 'g');
L(3) = plot(nan, nan, 'r');
legend(L, {'Sask Drive', '110 Street', '83 Avenue'},'Location','northeast', 'Autoupdate', 'off')
% fill([125;125;200;200],[-8;12;12;-8],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
% fill([200;200;300;300],[-8;12;12;-8],'w','FaceAlpha',0.1,'EdgeAlpha', '1', 'Linestyle', ':');
fill([118;118;218;218],[-8;12;12;-8],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
fill([218;218;318;318],[-8;12;12;-8],'w','FaceAlpha',0.1,'EdgeAlpha', '1', 'Linestyle', ':');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              TARGET N1 ERP VARIABLES FOR JASP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ttest for targets
electrode = 15; %adjust the electrode number accordingly btn Fz and Pz
time_window = find(EEG.times>118,1)-1:find(EEG.times>218,1)-2;
%[h p ci stat] = ttest(squeeze(mean(erp_out(time_window,1,electrode,1,:),1)),squeeze(mean(erp_out(time_window,1,electrode,3,:),1)),.05,'both',1)

% TARGET VARIABLES FOR JASP
N1_targ_sask = squeeze(mean(erp_out(time_window,1,electrode,1,:),1));
N1_targ_110 = squeeze(mean(erp_out(time_window,1,electrode,2,:),1));
N1_targ_83 = squeeze(mean(erp_out(time_window,1,electrode,3,:),1));
N1_targ_conds = [N1_targ_sask, N1_targ_110, N1_targ_83];
% writematrix(N1_targ_conds, 'Fz_N1_targets.csv')
writematrix(N1_targ_conds, 'Pz_N1_targets.csv')
clear time_window;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              TARGET P2 ERP VARIABLES FOR JASP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


electrode = 15; %adjust the electrode number accordingly btn Fz and Pz
time_window = find(EEG.times>218,1)-1:find(EEG.times>318,1)-2;
P2_targ_sask = squeeze(mean(erp_out(time_window,1,electrode,1,:),1));
P2_targ_110 = squeeze(mean(erp_out(time_window,1,electrode,2,:),1));
P2_targ_83 = squeeze(mean(erp_out(time_window,1,electrode,3,:),1));
P2_targ_conds = [P2_targ_sask, P2_targ_110, P2_targ_83];
% writematrix(P2_targ_conds,'Fz_P2_targets.csv')
writematrix(P2_targ_conds,'Pz_P2_targets.csv')
clear time_window;

%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         STANDARDS ERP FIGURES + STATS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Standards on same axis
electrode = 15;
figure;
boundedline(EEG.times,squeeze(mean(erp_out(:,2,electrode,1,:),5)), squeeze(std(erp_out(:,2,electrode,1,:),[],5))./sqrt(nsubs),'b',...
    EEG.times,squeeze(mean(erp_out(:,2,electrode,2,:),5)), squeeze(std(erp_out(:,2,electrode,2,:),[],5))./sqrt(nsubs),'g', ...
    EEG.times,squeeze(mean(erp_out(:,2,electrode,3,:),5)), squeeze(std(erp_out(:,2,electrode,3,:),[],5))./sqrt(nsubs),'r');
set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Standard Tones');
xlabel('Time (ms)');
ylabel('Voltage (uV)');

L(1) = plot(nan, nan, 'b');
L(2) = plot(nan, nan, 'g');
L(3) = plot(nan, nan, 'r');
legend(L, {'Sask Drive', '110 Street', '83 Avenue'},'Location','northeast', 'Autoupdate', 'off')
% fill([125;125;200;200],[-8;12;12;-8],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
% fill([200;200;300;300],[-8;12;12;-8],'w','FaceAlpha',0.1,'EdgeAlpha', '1', 'Linestyle', ':');
fill([118;118;218;218],[-8;12;12;-8],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
fill([218;218;318;318],[-8;12;12;-8],'w','FaceAlpha',0.1,'EdgeAlpha', '1', 'Linestyle', ':');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              STANDARD N1 ERP VARIABLES FOR JASP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
electrode = 15;
time_window = find(EEG.times>118,1)-1:find(EEG.times>218,1)-2;
N1_stan_sask = squeeze(mean(erp_out(time_window,2,electrode,1,:),1));
N1_stan_110 = squeeze(mean(erp_out(time_window,2,electrode,2,:),1)); 
N1_stan_83 = squeeze(mean(erp_out(time_window,2,electrode,3,:),1));
N1_stan_conds = [N1_stan_sask, N1_stan_110, N1_stan_83];
%writematrix(N1_stan_conds, 'Fz_N1_standards.csv')
 writematrix(N1_stan_conds, 'Pz_N1_standards.csv')
clear time_window;
clear electrode;
%%
% time_window = find(EEG.times>225,1)-1:find(EEG.times>325,1)-2;
%[h p ci stat] = ttest(squeeze(mean(erp_out(time_window,2,electrode,3,:),1)),squeeze(mean(erp_out(time_window,2,electrode,1,:),1)),.05,'both',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              STANDARD P2 ERP VARIABLES FOR JASP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
electrode = 15;
time_window = find(EEG.times>218,1)-1:find(EEG.times>318,1)-2;
P2_stan_sask = squeeze(mean(erp_out(time_window,2,electrode,1,:),1));
P2_stan_110 = squeeze(mean(erp_out(time_window,2,electrode,2,:),1));
P2_stan_83 = squeeze(mean(erp_out(time_window,2,electrode,3,:),1));
P2_stan_conds = [P2_stan_sask, P2_stan_110, P2_stan_83];
% writematrix(P2_stan_conds, 'Fz_P2_standards.csv')
writematrix(P2_stan_conds, 'Pz_P2_standards.csv')
clear time_window;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%                MANUSCRIPT TOPOGRAPHIES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%difference topographies for P3 figure
erp_diff_out = squeeze(erp_out(:,1,:,:,:)-erp_out(:,2,:,:,:));
time_window = find(EEG.times>355,1)-1:find(EEG.times>505,1)-2;
figure('Color',[1 1 1]);
for i_cond = 1:nconds
    subplot(1,nconds,i_cond);
    set(gca,'Color',[1 1 1]);
    temp = mean(mean(erp_diff_out(time_window,:,i_cond,:),4),1);
    temp(16:18) = NaN;
    topoplot(temp,'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4])
    title(conds_lab{i_cond});
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');
    
end
%%
%Target & Standards topographies
time_window = find(EEG.times>218,1)-1:find(EEG.times>318,1)-2;
figure('Color',[1 1 1]);
subplot(1,3,1);
set(gca,'Color',[1 1 1]);
temp1 = mean(mean(erp_out(time_window,1,:,1,:),5),1);
temp1(16:18) = NaN;
topoplot(temp1,'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4])
title('Sask Drive');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Targets');

subplot(1,3,2);
set(gca,'Color',[1 1 1]);
temp2 = mean(mean(erp_out(time_window,1,:,2,:),5),1);
temp2(16:18) = NaN;
topoplot(temp2,'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4])
title('110 Street');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Targets');

subplot(1,3,3);
set(gca,'Color',[1 1 1]);
temp3 = mean(mean(erp_out(time_window,1,:,3,:),5),1);
temp3(16:18) = NaN;
topoplot(temp3,'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4])
title('83 Avenue');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Targets');

%standards
figure('Color',[1 1 1]);
subplot(1,3,1);
set(gca,'Color',[1 1 1]);
temp4 = mean(mean(erp_out(time_window,2,:,1,:),5),1);
temp4(16:18) = NaN;
topoplot(temp4,'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4])
title('Sask Drive');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Standards');

subplot(1,3,2);
set(gca,'Color',[1 1 1]);
temp5 = mean(mean(erp_out(time_window,2,:,2,:),5),1);
temp5(16:18) = NaN;
topoplot(temp5,'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4])
title('110 Street');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Standards');

subplot(1,3,3);
set(gca,'Color',[1 1 1]);
temp6 = mean(mean(erp_out(time_window,2,:,3,:),5),1);
temp6(16:18) = NaN;
topoplot(temp6,'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4])
title('83 Avenue');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Standards');

clear time_window 



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% TIME WINDOW + PEAK AVERAGES TO DETERMINE STATS WINDOW.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P3 PEAK (DIFFERENCE-WAVES)
%trying out on difference
erp_diff_out = squeeze(erp_out(:,1,:,:,:)-erp_out(:,2,:,:,:));
electrode = 15;
figure;
boundedline(EEG.times,squeeze(mean(mean(erp_diff_out(:,electrode,:,:),4),3)), squeeze(std(0))./sqrt(nsubs),'m'),...
    set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Difference-Wave Grand ERPs');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
fill([355;355;505;505],[-8;12;12;-8],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
clear electrode;

% MMN/N2b % NOT INCLUDING IN MANUSCRIPT UNLESS SIGNIFICANT
% electrode = 13;
% figure;
% boundedline(EEG.times,squeeze(mean(mean(erp_diff_out(:,electrode,:,:),4),3)), squeeze(std(0))./sqrt(nsubs),'m'),...
% set(gca,'Color',[1 1 1]);
% set(gca,'YDir','reverse');
% axis tight; ylim([-12 16]);
% line([-200 1000],[0 0],'color','k');
% line([0 0],[-2.5 8],'color','k');
% title('Difference-Wave Grand ERPs');
% xlabel('Time (ms)');
% ylabel('Voltage (uV)');
% clear electrode;

% N1 + P2 PEAK/TIME WINDOW - TARGETS
electrode = 13;
figure;
boundedline(EEG.times,squeeze(mean(mean(erp_out(:,1,electrode,:,:),4),5)), squeeze(std(0))./sqrt(nsubs),'m'),...
    set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Target Grand ERPs');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
 fill([118;118;218;218],[-8;12;12;-8],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
 fill([218;218;318;318],[-8;12;12;-8],'w','FaceAlpha',0.1,'EdgeAlpha', '1', 'Linestyle', ':');
clear electrode;

% N1 + P2 PEAK/TIME WINDOW - STANDARDS
electrode = 13;
figure;
boundedline(EEG.times,squeeze(mean(mean(erp_out(:,2,electrode,:,:),4),5)), squeeze(std(0))./sqrt(nsubs),'m'),...
    set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Standard Grand ERPs');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
 fill([118;118;218;218],[-8;12;12;-8],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
 fill([218;218;318;318],[-8;12;12;-8],'w','FaceAlpha',0.1,'EdgeAlpha', '1', 'Linestyle', ':');
clear electrode;

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                  BARGRAPHS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% N1  FZ standards and targets

%%%%%%%%%%%%%%%%%%%%%%%%%%
electrode = 13;
time_window = find(EEG.times>118,1)-1:find(EEG.times>218,1)-2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%STANDARDS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hmean = mean(squeeze(mean(erp_out(time_window,2,electrode,1,:),5)));
imean = mean(squeeze(mean(erp_out(time_window,2,electrode,2,:),5)));
lmean = mean(squeeze(mean(erp_out(time_window,2,electrode,3,:),5)));

hsd = mean(squeeze(std(erp_out(time_window,2,electrode,1,:),[],5))./sqrt(nsubs));
isd = mean (squeeze(std(erp_out(time_window,2,electrode,2,:),[],5))./sqrt(nsubs));
lsd = mean(squeeze(std(erp_out(time_window,2,electrode,3,:),[],5))./sqrt(nsubs));

n1_mean_fz = [hmean; imean; lmean]';
n1_sd_fz = [hsd; isd; lsd]';

conds_plot = {'Heavy'; 'Intermediate';'Low'}; 
figure;
subplot (2,4,1)
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(n1_mean_fz,n1_sd_fz);
ylim([-6 0])
%ylabel('Standards Fz')
% xlabel('Condition')
title('Standards Fz')
legend(conds_plot)

clear hmean imean lmean hsd isd lsd n1_mean_fz n1_sd_fz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% TARGET N1 TIME WINDOW FZ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hmean = mean(squeeze(mean(erp_out(time_window,1,electrode,1,:),5)));
imean = mean(squeeze(mean(erp_out(time_window,1,electrode,2,:),5)));
lmean = mean(squeeze(mean(erp_out(time_window,1,electrode,3,:),5)));

hsd = mean(squeeze(std(erp_out(time_window,1,electrode,1,:),[],5))./sqrt(nsubs));
isd = mean (squeeze(std(erp_out(time_window,1,electrode,2,:),[],5))./sqrt(nsubs));
lsd = mean(squeeze(std(erp_out(time_window,1,electrode,3,:),[],5))./sqrt(nsubs));

tn1_mean_fz = [hmean; imean; lmean]';
tn1_sd_fz = [hsd; isd; lsd]';

conds_plot = {'Heavy'; 'Intermediate';'Low'}; 
subplot (2,4,2)
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(tn1_mean_fz,tn1_sd_fz);
ylim([-6 0])
%ylabel('Targets Fz')
% xlabel('Condition')
title('Targets Fz')
%legend(conds_plot)

clear hmean imean lmean hsd isd lsd n1_mean_fz n1_sd_fz electrode time_window


%%%%%%%%%%%%%%%%%%%%%
%
%
%        N1 PZ
%
%
%%%%%%%%%%%%%%%%%%%%%
electrode = 15;
time_window = find(EEG.times>118,1)-1:find(EEG.times>218,1)-2;

%STANDARDS
hmean = mean(squeeze(mean(erp_out(time_window,2,electrode,1,:),5)));
imean = mean(squeeze(mean(erp_out(time_window,2,electrode,2,:),5)));
lmean = mean(squeeze(mean(erp_out(time_window,2,electrode,3,:),5)));

hsd = mean(squeeze(std(erp_out(time_window,2,electrode,1,:),[],5))./sqrt(nsubs));
isd = mean (squeeze(std(erp_out(time_window,2,electrode,2,:),[],5))./sqrt(nsubs));
lsd = mean(squeeze(std(erp_out(time_window,2,electrode,3,:),[],5))./sqrt(nsubs));

n1_mean_pz = [hmean; imean; lmean]';
n1_sd_pz = [hsd; isd; lsd]';

conds_plot = {'Heavy'; 'Intermediate';'Low'}; 
subplot (2,4,3)
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(n1_mean_pz,n1_sd_pz);
ylim([-6 0])
%ylabel('Standards Pz')
% xlabel('Condition')
title('Standards Pz')
%legend(conds_plot)

clear hmean imean lmean hsd isd lsd n1_mean_fz n1_sd_fz

%TARGET N1 TIME WINDOW FZ

%Targets

hmean = mean(squeeze(mean(erp_out(time_window,1,electrode,1,:),5)));
imean = mean(squeeze(mean(erp_out(time_window,1,electrode,2,:),5)));
lmean = mean(squeeze(mean(erp_out(time_window,1,electrode,3,:),5)));

hsd = mean(squeeze(std(erp_out(time_window,1,electrode,1,:),[],5))./sqrt(nsubs));
isd = mean (squeeze(std(erp_out(time_window,1,electrode,2,:),[],5))./sqrt(nsubs));
lsd = mean(squeeze(std(erp_out(time_window,1,electrode,3,:),[],5))./sqrt(nsubs));

tn1_mean_pz = [hmean; imean; lmean]';
tn1_sd_pz = [hsd; isd; lsd]';

conds_plot = {'Heavy'; 'Intermediate';'Low'}; 
subplot (2,4,4)
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(tn1_mean_pz,tn1_sd_pz);
ylim([-6 0])
%ylabel('Targets Pz')
% xlabel('Condition')
title('Targets Pz')
%legend(conds_plot)

clear hmean imean lmean hsd isd lsd n1_mean_fz n1_sd_fz electrode time_window



%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    P2 FZ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_window = find(EEG.times>218,1)-1:find(EEG.times>318,1)-2;
electrode = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%STANDARDS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hmean = mean(squeeze(mean(erp_out(time_window,2,electrode,1,:),5)));
imean = mean(squeeze(mean(erp_out(time_window,2,electrode,2,:),5)));
lmean = mean(squeeze(mean(erp_out(time_window,2,electrode,3,:),5)));

hsd = mean(squeeze(std(erp_out(time_window,2,electrode,1,:),[],5))./sqrt(nsubs));
isd = mean (squeeze(std(erp_out(time_window,2,electrode,2,:),[],5))./sqrt(nsubs));
lsd = mean(squeeze(std(erp_out(time_window,2,electrode,3,:),[],5))./sqrt(nsubs));

P2_mean_fz = [hmean; imean; lmean]';
P2_sd_fz = [hsd; isd; lsd]';

conds_plot = {'Heavy'; 'Intermediate';'Low'}; 
subplot (2,4,5)
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(P2_mean_fz,P2_sd_fz);
ylim([-2 0])
%ylabel('Standards Pz')
% xlabel('Condition')
title('Standards Pz')
legend(conds_plot)

clear hmean imean lmean hsd isd lsd n1_mean_fz n1_sd_fz 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% TARGET 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hmean = mean(squeeze(mean(erp_out(time_window,1,electrode,1,:),5)));
imean = mean(squeeze(mean(erp_out(time_window,1,electrode,2,:),5)));
lmean = mean(squeeze(mean(erp_out(time_window,1,electrode,3,:),5)));

hsd = mean(squeeze(std(erp_out(time_window,1,electrode,1,:),[],5))./sqrt(nsubs));
isd = mean (squeeze(std(erp_out(time_window,1,electrode,2,:),[],5))./sqrt(nsubs));
lsd = mean(squeeze(std(erp_out(time_window,1,electrode,3,:),[],5))./sqrt(nsubs));

tp2_mean_fz = [hmean; imean; lmean]';
tp2_sd_fz = [hsd; isd; lsd]';

conds_plot = {'Heavy'; 'Intermediate';'Low'}; 
subplot (2,4,6)
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(tp2_mean_fz,tp2_sd_fz);
ylim([-5 0])
%ylabel('Targets Fz')
% xlabel('Condition')
title('Targets Fz')
%legend(conds_plot)

clear hmean imean lmean hsd isd lsd n1_mean_fz n1_sd_fz electrode



%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    P2 PZ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

electrode = 15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%STANDARDS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hmean = mean(squeeze(mean(erp_out(time_window,2,electrode,1,:),5)));
imean = mean(squeeze(mean(erp_out(time_window,2,electrode,2,:),5)));
lmean = mean(squeeze(mean(erp_out(time_window,2,electrode,3,:),5)));

hsd = mean(squeeze(std(erp_out(time_window,2,electrode,1,:),[],5))./sqrt(nsubs));
isd = mean (squeeze(std(erp_out(time_window,2,electrode,2,:),[],5))./sqrt(nsubs));
lsd = mean(squeeze(std(erp_out(time_window,2,electrode,3,:),[],5))./sqrt(nsubs));

p2_mean_pz = [hmean; imean; lmean]';
p2_sd_pz = [hsd; isd; lsd]';

conds_plot = {'Heavy'; 'Intermediate';'Low'}; 
subplot (2,4,7)
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(p2_mean_pz,p2_sd_pz);
ylim([0 2])
%ylabel('Standards Pz')
% xlabel('Condition')
title('Standards Pz')
%legend(conds_plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% TARGET 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hmean = mean(squeeze(mean(erp_out(time_window,1,electrode,1,:),5)));
imean = mean(squeeze(mean(erp_out(time_window,1,electrode,2,:),5)));
lmean = mean(squeeze(mean(erp_out(time_window,1,electrode,3,:),5)));

hsd = mean(squeeze(std(erp_out(time_window,1,electrode,1,:),[],5))./sqrt(nsubs));
isd = mean (squeeze(std(erp_out(time_window,1,electrode,2,:),[],5))./sqrt(nsubs));
lsd = mean(squeeze(std(erp_out(time_window,1,electrode,3,:),[],5))./sqrt(nsubs));

tp2_mean_pz = [hmean; imean; lmean]';
tp2_sd_pz = [hsd; isd; lsd]';

conds_plot = {'Heavy'; 'Intermediate';'Low'}; 
subplot (2,4,8)
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(tp2_mean_pz,tp2_sd_pz);
ylim([-1.4 .9])
%ylabel('Targets Pz')
% xlabel('Condition')
title('Targets Pz')
%legend(conds_plot)

clear hmean imean lmean hsd isd lsd n1_mean_fz n1_sd_fz electrode


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXPLORATORY ANALYSES 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%individual subs

for i_sub = 1:nsubs

         figure
        boundedline(EEG.times,squeeze(erp_out(:,1,electrode,:,i_sub)),squeeze(std(0))./sqrt(nsubs),'r',...
            EEG.times,squeeze(erp_out(:,2,electrode,:,i_sub)),squeeze(std(0))./sqrt(nsubs),'k');
        set(gca,'Color',[1 1 1]);
        set(gca,'YDir','reverse');
        axis tight; ylim([-25 25]);
        line([-200 1000],[0 0],'color','k');
        line([0 0],[-2.5 8],'color','k');
        title(subs{i_sub});
        xlabel('Time (ms)');
        ylabel('Voltage (uV)');
 end

%%
%%
% %SINGLE ERPS FOR INDIVIDUAL COPYING
% 
% for i_cond = 1:nconds
%     switch i_cond
%         case 1
%             colour = 'b';
%         case 2
%             colour = 'g';
%         case 3
%             colour = 'r';
%     end
%     
%     figure;
%     boundedline(EEG.times,squeeze(mean(erp_out(:,1,electrode,i_cond,:),5)),squeeze(std(erp_out(:,1,electrode,i_cond,:),[],5))./sqrt(nsubs),colour,...
%         EEG.times,squeeze(mean(erp_out(:,2,electrode,i_cond,:),5)),squeeze(std(erp_out(:,2,electrode,i_cond,:),[],5))./sqrt(nsubs),'k');
%     set(gca,'Color',[1 1 1]);
%     set(gca,'YDir','reverse');
%     if i_cond == 2 || 3
%         legend('Targets','Standards','Location','NorthEast');
%     elseif i_cond == 1
%         legend('Targets','Standards','Location','NorthEast');
%     end
%     axis tight; ylim([-8 12]);
%     line([-200 1000],[0 0],'color','k');
%     line([0 0],[-2.5 8],'color','k');
%     title(conds{i_cond});
%     xlabel('Time (ms)');
%     ylabel('Voltage (uV)');
% end

%%
%%

% for i_set = 1:48; trial_count(i_set) = ALLEEG(i_set).trials; end
% trial_count = reshape(trial_count,[2,3,8]);
% min(trial_count,[],3)
% mean(trial_count,3)
% max(trial_count,[],3)
%
% %mean and sd
% mean(mean(erp_diff_out(time_window,7,1:3,:),1),4)
% std(mean(erp_diff_out(time_window,7,1:3,:),1),[],4)
%
%
%
% ttest of each condition
% [h p ci stat] = ttest(squeeze(mean(erp_diff_out(time_window,7,1,:),1)),0,.05,'right',1)
% [h p ci stat] = ttest(squeeze(mean(erp_diff_out(time_window,7,2,:),1)),0,.05,'right',1)
% [h p ci stat] = ttest(squeeze(mean(erp_diff_out(time_window,7,3,:),1)),0,.05,'right',1)
% [h p ci stat] = ttest(squeeze(mean(erp_diff_out(time_window,7,4,:),1)),0,.05,'right',1)
% 
% 
% 
% eeglab redraw

