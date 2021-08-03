clear all
close all
clc

addpath('StructFiles/h_ranges')

%% 1st range
load('Synth-start-1-h-10-pm-30-N-90.mat')
T1 = INFO.TABLE;
N1 = INFO.Training_set;

load('Synth-start-1-h-10-pm-30-N-100.mat')
T2 = INFO.TABLE;
N2 = INFO.Training_set;

load('Synth-start-1-h-10-pm-30-N-150.mat')
T3 = INFO.TABLE;
N3 = INFO.Training_set;

load('Synth-start-1-h-10-pm-30-N-200.mat')
T4 = INFO.TABLE;
N4 = INFO.Training_set;

load('Synth-start-1-h-10-pm-30-N-300.mat')
T5 = INFO.TABLE;
N5 = INFO.Training_set;

load('Synth-start-1-h-10-pm-30-N-400.mat')
T6 = INFO.TABLE;
N6 = INFO.Training_set;

load('Synth-start-1-h-10-pm-30-N-600.mat')
T7 = INFO.TABLE;
N7 = INFO.Training_set;

load('Synth-start-1-h-10-pm-30-N-800.mat')
T8 = INFO.TABLE;
N8 = INFO.Training_set;


%% 2nd range

% load('Synth-start-6-h-10-pm-30-N-90.mat')
% T1 = INFO.TABLE;
% N1 = INFO.Training_set;
% 
% load('Synth-start-6-h-10-pm-30-N-100.mat')
% T2 = INFO.TABLE;
% N2 = INFO.Training_set;
% 
% load('Synth-start-6-h-10-pm-30-N-150.mat')
% T3 = INFO.TABLE;
% N3 = INFO.Training_set;
% 
% load('Synth-start-6-h-10-pm-30-N-200.mat')
% T4 = INFO.TABLE;
% N4 = INFO.Training_set;
% 
% load('Synth-start-6-h-10-pm-30-N-300.mat')
% T5 = INFO.TABLE;
% N5 = INFO.Training_set;
% 
% load('Synth-start-6-h-10-pm-30-N-400.mat')
% T6 = INFO.TABLE;
% N6 = INFO.Training_set;
% 
% load('Synth-start-6-h-10-pm-30-N-600.mat')
% T7 = INFO.TABLE;
% N7 = INFO.Training_set;
% 
% load('Synth-start-6-h-10-pm-30-N-800.mat')
% T8 = INFO.TABLE;
% N8 = INFO.Training_set;

% load('Info-WSA-250-27-20-10.mat')
% T9 = INFO.TABLE;
% N9 = INFO.Training_set;
% 
% load('Info-WSA-300-27-20-10.mat')
% T10 = INFO.TABLE;
% N10 = INFO.Training_set;
% 
% load('Info-WSA-400-27-20-10.mat')
% T11 = INFO.TABLE;
% N11 = INFO.Training_set;
% 
% load('Info-WSA-500-27-20-10.mat')
% T12 = INFO.TABLE;
% N12 = INFO.Training_set;
% 
% load('Info-WSA-600-27-20-10.mat')
% T13 = INFO.TABLE;
% N13 = INFO.Training_set;
% 
% load('Info-WSA-800-27-20-10.mat')
% T14 = INFO.TABLE;
% N14 = INFO.Training_set;
% 
% load('Info-WSA-1000-27-20-10.mat')
% T15 = INFO.TABLE;
% N15 = INFO.Training_set;


Nvec = [N1 N2 N3 N4 N5 N6 N7 N8];

p1 = table2array(T1(:,5));
p2 = table2array(T2(:,5));
p3 = table2array(T3(:,5));
p4 = table2array(T4(:,5));
p5 = table2array(T5(:,5));
p6 = table2array(T6(:,5));
p7 = table2array(T7(:,5));
p8 = table2array(T8(:,5));
% p9 = table2array(T9(:,5));
% p10 = table2array(T10(:,5));
% p11 = table2array(T11(:,5));
% p12 = table2array(T12(:,5));
% p13 = table2array(T13(:,5));
% p14 = table2array(T14(:,5));
% p15 = table2array(T15(:,5));

e1 = table2array(T1(:,2));
e2 = table2array(T2(:,2));
e3 = table2array(T3(:,2));
e4 = table2array(T4(:,2));
e5 = table2array(T5(:,2));
e6 = table2array(T6(:,2));
e7 = table2array(T7(:,2));
e8 = table2array(T8(:,2));
% e9 = table2array(T9(:,2));
% e10 = table2array(T10(:,2));
% e11 = table2array(T11(:,2));
% e12 = table2array(T12(:,2));
% e13 = table2array(T13(:,2));
% e14 = table2array(T14(:,2));
% e15 = table2array(T15(:,2));

pvec = [p1 p2 p3 p4 p5 p6 p7 p8];
evec = [e1 e2 e3 e4 e5 e6 e7 e8];

%% PLOTS
figure
plot(Nvec, pvec','LineWidth',1.5)
legend('BVIC(r=2)','BVIC(r=8)','AIC','BIC','AICc')
xlabel('N')
ylabel('Average order selected')
title('Orders selected by criteria')
figure
plot(Nvec, evec,'LineWidth',1)
xlabel('N')
ylabel('Average MSFE')
title('average MSFE for each sample size')
legend('BVIC(r=2)','BVIC(r=8)','AIC','BIC','AICc')
