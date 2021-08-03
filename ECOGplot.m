clear all
close all
clc

addpath('Seizure_DATA')
load('clinical_metadata.mat')
load('HUP68-ictal-block-1.mat')
range = 12;
X = evData(range,:);
c = 100+(range-range(1)+1)*1000;
t = linspace(0,round(length(X)/500),length(X));
Xp = X+c';

channels = subject(3).Channels(range);


% wo = 60/(500/2);
% bw = wo/35;
% Ab = 3;
% [fnum,fden] = iirnotch(wo,bw,Ab);
% Y = filter(fnum,fden,X);


figure
plot(Xp','k')

% xlim([10 20])
set(gca,'ytick',c','FontSize',11)
% set(gca, 'xticklabels',[])
yticklabels(channels)
% xticklabels({'t_0', 't_0 + 1','t_0 + 2','t_0 + 3','t_0 + 4','t_0 + 5'})
xlabel('Time','FontSize',18)
ylabel('Channels','FontSize',18)
hold on
line([14.5 15.5],[-1 -1],'Color','k','LineWidth',2)
grid on
% text(14.75,-400,'1s','FontSize',14,'FontWeight','bold')



% figure
% pspectrum(X)
% hold on
% pspectrum(Y)