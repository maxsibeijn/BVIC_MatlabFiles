clear all
close all
clc

addpath('StructFiles')
% load('Info-ECoG-35-1-50-200.mat')
% T1 = INFO.TABLE;
% h1 = INFO.Horizon;
% 
% load('Info-ECoG-35-2-50-200.mat')
% T2 = INFO.TABLE;
% h2 = INFO.Horizon;
% 
% load('Info-ECoG-35-4-50-200.mat')
% T3 = INFO.TABLE;
% h3 = INFO.Horizon;
% 
% load('Info-ECoG-35-6-50-200.mat')
% T4 = INFO.TABLE;
% h4 = INFO.Horizon;
% 
% load('Info-ECoG-35-10-50-200.mat')
% T5 = INFO.TABLE;
% h5 = INFO.Horizon;
% 
% load('Info-ECoG-35-14-50-200.mat')
% T6 = INFO.TABLE;
% h6 = INFO.Horizon;
% 
% load('Info-ECoG-35-18-50-200.mat')
% T7 = INFO.TABLE;
% h7 = INFO.Horizon;
% 
% load('Info-ECoG-35-22-50-200.mat')
% T8 = INFO.TABLE;
% h8 = INFO.Horizon;
% load('Info-ECoG-35-26-50-200.mat')
% T9 = INFO.TABLE;
% h9 = INFO.Horizon;
% 
% load('Info-ECoG-35-30-50-200.mat')
% T10 = INFO.TABLE;
% h10 = INFO.Horizon;

%% single point forecast

load('ECoG-start-1-h-25-pm-70-N-220.mat')
T1 = INFO.TABLE;
h1 = 1;

load('ECoG-start-2-h-25-pm-70-N-220.mat')
T2 = INFO.TABLE;
h2 = 2;

load('ECoG-start-4-h-25-pm-70-N-220.mat')
T3 = INFO.TABLE;
h3 = 4;

load('ECoG-start-6-h-25-pm-70-N-220.mat')
T4 = INFO.TABLE;
h4 = 6;

load('ECoG-start-8-h-25-pm-70-N-220.mat')
T5 = INFO.TABLE;
h5 = 8;

load('ECoG-start-10-h-25-pm-70-N-220.mat')
T6 = INFO.TABLE;
h6 = 10;

load('ECoG-start-12-h-25-pm-70-N-220.mat')
T7 = INFO.TABLE;
h7 = 12;

load('ECoG-start-14-h-25-pm-70-N-220.mat')
T8 = INFO.TABLE;
h8 = 14;
load('ECoG-start-16-h-25-pm-70-N-220.mat')
T9 = INFO.TABLE;
h9 = 16;

load('ECoG-start-18-h-25-pm-70-N-220.mat')
T10 = INFO.TABLE;
h10 = 18;

load('ECoG-start-20-h-25-pm-70-N-220.mat')
T11 = INFO.TABLE;
h11 = 20;

load('ECoG-start-22-h-25-pm-70-N-220.mat')
T12 = INFO.TABLE;
h12 = 22;

load('ECoG-start-24-h-25-pm-70-N-220.mat')
T13 = INFO.TABLE;
h13 = 24;


Nvec = [h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 h11 h12 h13];

p1 = table2array(T1(:,5));
p2 = table2array(T2(:,5));
p3 = table2array(T3(:,5));
p4 = table2array(T4(:,5));
p5 = table2array(T5(:,5));
p6 = table2array(T6(:,5));
p7 = table2array(T7(:,5));
p8 = table2array(T8(:,5));
p9 = table2array(T9(:,5));
p10 = table2array(T10(:,5));
p11 = table2array(T11(:,5));
p12 = table2array(T12(:,5));
p13 = table2array(T13(:,5));

e1 = table2array(T1(:,2));
e2 = table2array(T2(:,2));
e3 = table2array(T3(:,2));
e4 = table2array(T4(:,2));
e5 = table2array(T5(:,2));
e6 = table2array(T6(:,2));
e7 = table2array(T7(:,2));
e8 = table2array(T8(:,2));
e9 = table2array(T9(:,2));
e10 = table2array(T10(:,2));
e11 = table2array(T11(:,2));
e12 = table2array(T12(:,2));
e13 = table2array(T13(:,2));

idx = 3;
pvec = [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13];
% evec = [e1/e1(idx) e2/e2(idx) e3/e3(idx) e4/e4(idx) e5/e5(idx) e6/e6(idx) e7/e7(idx) e8/e8(idx) e9/e9(idx) e10/e10(idx) e11/e11(idx) e12/e12(idx) e13/e13(idx)];
 evec = [e1 e2 e3 e4 e5 e6 e7 e8 e9 e10 e11 e12 e13];

% evec(4,:) = [];
% pvec(4,:) = [];
%% PLOTS
figure
plot(Nvec, pvec','LineWidth',1.5)
legend('BVIC(r=2)','BVIC(r=8)','AIC','BIC','AICc')
xlabel('h')
ylabel('Average order selected')
title('Orders selected by criteria')
figure
plot(Nvec, evec','LineWidth',1.5)
xlabel('h')
ylabel('Average MSFE')
title('Average MSFE for each prediction horizon ')
legend('BVIC(r=2)','BVIC(r=8)','AIC','BIC','AICc')