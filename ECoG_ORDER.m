clear all 
close all
clc

% 500hz reasonable assumption?
load('HUP68-ictal-block-1.mat')

XX = evData(1,:)';
 X = evData(1,13001:14000)';
% X = evData(1,end-2000:end)';
% 
p = 48;
np = 4:4:p;
qq = numel(np);
N = length(X);

tic
for i  = 1:qq
    
Md(i) = arima(np(i),0,0);
[~,~,logL(i),~] = estimate(Md(i),X,'Display','off');


end
toc

[aic, bic] = aicbic(logL,np,N);

[~,pa] = min(aic);
[~,pb] = min(bic);

pa = np(pa)
pb = np(pb)
