clear all
close all
clc

addpath('Seizure_DATA')

rng default

DataType = 1; % 0 if generated(small order), 1 if generated(big order), 2 if real data

    
order =4; %Order of the generated model
mm = 5; % number of steps used for forecasting

if DataType == 0
    N = 100; % size of window
    MonteCarloLoops = 300;
    Windows = 100; % number of windows
elseif DataType == 1
    N = 200; % size of window
    MonteCarloLoops = 10;
    Windows = 100; % number of windows
elseif DataType == 2
    load('HUP68-interictal-block-1.mat') % load the specific dataset to be analysed
    Ndata = 10000; %floor(length(evData)/1000)*1000; % Usable length of Ecog timeseries
    Tstart = 20000;
    N = 200; % size of window
    MonteCarloLoops = size(evData,1);
    Windows = Ndata/N; % number of windows
end


tic
counter = 0;

%% Variable PreAllocation
MSE_B = zeros(Windows,MonteCarloLoops);
MSE_F = zeros(Windows,MonteCarloLoops);
MSE_AIC = zeros(Windows,MonteCarloLoops);
MSE_ML = zeros(Windows,MonteCarloLoops);
V_B = zeros(Windows,MonteCarloLoops);
V_F = zeros(Windows,MonteCarloLoops);
V_AIC = zeros(Windows,MonteCarloLoops);
V_ML = zeros(Windows,MonteCarloLoops);
p1 = zeros(Windows,MonteCarloLoops);
p2 = zeros(Windows,MonteCarloLoops);
p3 = zeros(Windows,MonteCarloLoops);
MaxB = zeros(MonteCarloLoops,1);
MaxF = zeros(MonteCarloLoops,1);
MaxAIC = zeros(MonteCarloLoops,1);
MaxML = zeros(MonteCarloLoops,1);
MinB = zeros(MonteCarloLoops,1);
MinF = zeros(MonteCarloLoops,1);
MinAIC = zeros(MonteCarloLoops,1);
MinML = zeros(MonteCarloLoops,1);


%% START OF LOOP
for mc = 1:MonteCarloLoops
    

if DataType < 1.5
[X,XM,phi] = TSgen(N,order, Windows);

else
X = evData(mc,1+Tstart:Ndata+Tstart);
X = (X-mean(X))/std(X); %standardization
% 
XM = zeros(N,Ndata/N);
for ii = 1:Ndata/N
    XM(:,ii) = X((ii-1)*N+1:ii*N); 
end
end

%% Model
for zz = 1:Windows
    XL = XM(:,zz); % The 100 samples belonging to the zz'th window
    Xloop =  (XL- mean(XL))/std(XL);
    NT = 0.5*N+round(2*(order)); % Determines the size of forecast set 
    XT = zeros(2*NT-N,1);
    XT = Xloop(N-NT+1:NT);

    XT2 = Xloop(N-NT-mm+1:NT);
    
    
    
% function [mlogL, mmse, mvar, acf] = bcast(Xloop,XT,N,NT,mm,order)


p = 1:order+5;%ceil(.2*order);
idx = 1:numel(p);
% acf = autocorr(XT,'NumLags',max(p)+mm);
acf2 = autocorr(XT2,'NumLags',max(p)+mm);
% logL = zeros(length(p),1);
biasX = zeros(length(p),mm);
varX = zeros(length(p),mm);
        

for i = 1:numel(p)
%     
%     Mdl(i) = arima(p(i),0,0);
%     [MD,~,logL3(i,zz),~] = estimate(Mdl(i),XT,'Display','off');

    % LogLikelihood computation by OLS: 
    %Choose Xloop for in-sample forecast/ Choose XT for out-of-sample
    %forecast
    Xk = XT2; % (XT-mean(XT))/std(XT); 
    rL = acf2(2:p(i)+1);
    R_n = toeplitz(acf2(1:p(i))); %AUTOCORRELATION MATRIX R
%      XLSc = toeplitz(Xk(p(i):end-1),flip(Xk(1:p(i))));
    XLS = hankel(flip(Xk(p(i):end-1)),flip(Xk(1:p(i)))); % tranposed hankel matrix
    YLS = flip(Xk(p(i)+1:end));
    phiLS = (XLS'*XLS)\XLS'*YLS;
    phi_L = R_n\rL;
 
% %     
% 
% 
%     
%     [logL(i), sigma2(i), ISL, logL2(i)] = loglikelihood(Xk,phiLS,p(i)); 

    ISL = length(Xk)-max(p);             % in-sample size
%     eps = zeros(ISL,1);
    for k = 1:ISL
    eps(k,i) = Xk(p(i)+k) - phiLS'*flip(Xk(k:p(i)+k-1));
    end
    sigma2(i) = (ISL)\sum(eps(:,i).^2);
%     logL(i,zz) = -ISL/2*(log(2*pi*sigma2(i)')-1);
    logL(i,zz) = -ISL*log(sigma2(i)');

%     logL2(i,zz) = -ISL/2*log(2*pi*sigma2(i)')-ISL/2;
%    


%     logL2(i,zz) = -ISL/2*(1/(var(Xk))*sigma2(i)');
   
    for m = 1:mm
        rho_n = acf2(m+1:m+p(i));        %AUTOCORRELATION FUNCTION rho
        phi_f = R_n\rho_n;              % PHI autoregressive coefficient vector
        errorX(i,m) = Xloop(N-NT-m+1) - phi_f'*XT(1:p(i));
        varX(i,m)   = var(XT)*(1- phi_f'*rho_n);
    end
    
    
    
    b2X(i,1) = norm(errorX(i,:))^2/mm;
    vX(i,1)  = sum(varX(i,:))/mm;
end

% 
% logL2 = logL2';

mlogL = logL(:,zz);
mmse = mean(b2X,2);
mvar = mean(vX,2);
   
% end
    
% [mlogL, mmse, mvar, acf] = bcast(Xloop,XT,N,NT,mm,order);
    



%% Criterion

[~,pml] = max(mlogL);
[~,ptest] = min(mlogL);

[~,pmse] = min(mmse);
[~,pvar] = min(mvar);

fvec = [zeros(size(mmse)), mmse/mmse(pml), mvar/mvar(pml)];

% BVIC parameters
Beta = 5;
Gamma = 1;


[obj, pmin] = bvic(fvec, Beta, Gamma,p); % the criterion
aic2 = -mlogL + 2*p';
[aic, bic] = aicbic(mlogL/2,p',ISL);
aicc = aic + (2*p'.^2+2*p')./(ISL - p' - 1);

[~,paic] = min(aic);
paic = p(paic);

[~,pc] = min(aicc);
pc = p(pc);
% 
% for a = 1:S+1
%     for b = 1:S+1
%         W = [1; 10/S*(b-1); 10/S*(a-1)];
%         [obj(a,b), pmin(a,b)] = min(fvec*W);
%         obj(a,b) = obj(a,b)/(10/S*(a+b-2)+1);
%         ax = linspace(0,10,S+1);
%     end
% end 



%% Forecast

pf = p(pmin);
pml = p(pml);


[f_err, f_aic, f_ml, f_c, Va, Vaic, Vml, Vc] = fcast(Xloop,XT,acf2,pf,paic,pml,pc,mm,NT);


%% Validation test
MSE_F(zz,mc) = norm(f_err)^2/mm;
MSE_B(zz,mc) = mmse(pmin);
MSE_AIC(zz,mc) = norm(f_aic)^2/mm;
MSE_ML(zz,mc) = norm(f_ml)^2/mm;
MSE_C(zz,mc) = norm(f_c)^2/mm;

V_F(zz,mc) = mean(Va); 
V_B(zz,mc) = mvar(pmin);
V_AIC(zz,mc) = mean(Vaic);
V_ML(zz,mc) = mean(Vml);
V_C(zz,mc) = mean(Vc);

CB = 1.96*sqrt(Va);


p1(zz,mc) = pf;
p2(zz,mc) = paic;
p3(zz,mc) = pml;
p4(zz,mc) = pc;
end


% 
% figure
% histogram(MSE_B,50)
% title('Backward MSE')
% figure
% histogram(MSE_F,50)
% title('Forward MSE of BVIC')
% figure
% histogram(MSE_AIC,50)
% title('Forward MSE of AIC')

% figure
% ecdf(MSE_B)
% hold on
% ecdf(MSE_F)
% ecdf(MSE_AIC)
% legend('Back','Forward','AIC')
MaxB(mc) = max(MSE_B(:,mc));
MaxF(mc) = max(MSE_F(:,mc));
MaxAIC(mc) = max(MSE_AIC(:,mc));
MaxML(mc) = max(MSE_ML(:,mc));
MinB(mc) = min(MSE_B(:,mc));
MinF(mc) = min(MSE_F(:,mc));
MinAIC(mc) = min(MSE_AIC(:,mc));
MinML(mc) = min(MSE_ML(:,mc));


[pt(mc),h(mc)] = signrank(MSE_F(:,mc),MSE_C(:,mc));



counter = counter+1;
if counter == ceil(.33*MonteCarloLoops)
    disp('one thirds done.')
end
if counter == ceil(.66*MonteCarloLoops)
   disp('two thirds done..') 
end
if counter == ceil(.90*MonteCarloLoops)
    disp('Almost done...')
end

end

%% Distributions

% One vector containng all MSE estimates
mseB = MSE_B(:);
mseF = MSE_F(:);
mseAIC = MSE_AIC(:);
mseML = MSE_ML(:);
mseC = MSE_C(:);

% vector containing all theoretical variance estimates
tvarB = V_B(:);
tvarF = V_F(:);
tvarAIC = V_AIC(:);
tvarML = V_ML(:);
tvarC = V_C(:);



M = max([mseF; mseB; mseAIC; mseML]);
MF = max(mseF);
range = linspace(0,1,length(mseF));
Fdi = fitdist(mseF/M,'beta');
Fpdf = betapdf(range,Fdi.a,Fdi.b);
Fcdf = betacdf(range,Fdi.a,Fdi.b);




MB = max(mseB);
Bdi = fitdist(mseB/M,'beta');
Bpdf = betapdf(range,Bdi.a,Bdi.b);
Bcdf = betacdf(range,Bdi.a,Bdi.b);


MAIC = max(mseAIC);
AICdi = fitdist(mseAIC/M,'beta');
AICpdf = betapdf(range,AICdi.a,AICdi.b);
AICcdf = betacdf(range,AICdi.a,AICdi.b);


MML = max(mseML);
MLdi = fitdist(mseML/M,'beta');
MLpdf = betapdf(range,MLdi.a,MLdi.b);
MLcdf = betacdf(range,MLdi.a,MLdi.b);



psaveF = p1(:);
psaveAIC = p2(:);
psaveML = p3(:);
psaveC = p4(:);
toc





%% figures
% close all
% 
% % 
% figure
% BW1 = .01;
% subplot(4,1,1),
% h1 = histogram(mseB,'BinWidth',BW1);
%  axis([0 1 0 4000])
% title('Backcast MSE of BVIC (global)')
% 
% subplot(4,1,2),
% h2 = histogram(mseF,'BinWidth',BW1);
%  axis([0 1 0 4000])
% title('Forecast MSE of BVIC (global)')
% % xlabel('MSE')
% % hold on
% % 
% % plot(range,Bpdf)
% % 
% % 
% subplot(4,1,3),
% h3 = histogram(mseAIC,'BinWidth',BW1);
%  axis([0 1 0 4000])
% title('Forecast MSE of AIC (global)')
% % xlabel('MSE')
% % hold on
% % plot(range,Fpdf)
% subplot(4,1,4),
% h4 = histogram(mseML,'BinWidth',BW1);
%  axis([0 1 0 4000])
% title('Forecast MSE of ML (global)')
% xlabel('MSE')
% % hold on
% % plot(range,AICpdf)
% % plot(range,MLpdf)
% % legend('Back','Forward','AIC','ML')
% 
% 
% %FACTORS
% f1 = max(h2.Values)/max(Fpdf);
% f2 = max(h3.Values)/max(AICpdf);
% f3 = max(h4.Values)/max(MLpdf);
% f4 = max(h1.Values)/max(Bpdf);
% 
% 
% 
% 
% 
% % title('Backward MSE (global)')
% % xlabel('MSE')
% % hold on
% % 
% % plot(range,Bpdf)
% 
% 
% figure
% % title('Forward MSE of BVIC (global)')
% xlabel('MSE')
% hold on
% plot(range*M,Bpdf*f4,':k','LineWidth',2)
% plot(range*M,Fpdf*f1,'b','LineWidth',2)
% 
% 
% % title('Forward MSE of AIC (global)')
% xlabel('MSE')
% hold on
% plot(range*M,AICpdf*f2,'-.r','LineWidth',2)
% plot(range*M,MLpdf*f3,'--m','LineWidth',2)
% 
% legend('BVIC(back)','BVIC(fore)','AIC','ML')
% xlim([0 10])
% 
% 
% 
% 
% 
% 
% 
% 
% %% Histogram for maximum MSE
% % 
% % hm1 = histogram(MaxB,'BinWidth',.1);
% % hB = hm1.Values;
% % 
% % hm2 = histogram(MaxF,'BinWidth',.1);
% % hF = hm2.Values;
% % 
% % hm3 = histogram(MaxAIC,'BinWidth',.1);
% % hA = hm3.Values;
% % 
% % hm4 = histogram(MaxML,'BinWidth',.1);
% % hML = hm4.Values;
% 
% 
% 
figure
BW2 = .001;
hold on
subplot(4,1,1), histogram(MaxB,'BinWidth',BW2)
% axis([0 60 0 15])

title('BVIC (backcast)')

subplot(4,1,2), histogram(MaxF,'BinWidth',BW2)
%  axis([0 60 0 15])
title('BVIC (forecast)')

subplot(4,1,3), histogram(MaxAIC,'BinWidth',BW2)
%  axis([0 60 0 15])


title('AIC')

subplot(4,1,4), histogram(MaxML,'BinWidth',BW2)
%  axis([0 60 0 15])
xlabel('MSE')
title('Maximum Likelihood')
% 
% 
% 
% figure
% BW3 = .0002;
% hold on
% subplot(4,1,1), histogram(MinB,'BinWidth',BW3)
% axis([0 .005 0 15])
% title('BVIC (backcast)')
% 
% 
% subplot(4,1,2), histogram(MinF,'BinWidth',BW3)
% axis([0 .005 0 15])
% title('BVIC (forecast)')
% 
% subplot(4,1,3), histogram(MinAIC,'BinWidth',BW3)
% axis([0 .005 0 15])
% 
% 
% title('AIC')
% 
% subplot(4,1,4), histogram(MinML,'BinWidth',BW3)
% axis([0 0.005 0 15])
% xlabel('MSE')
% title('Maximum Likelihood')


% figure
% histogram(MinF,50)
% 
% histogram(MinAIC,50)
% 
% histogram(MinML,50)
% figure
% hold on
% ecdf(mseB)
% 
% ecdf(mseF)
% ecdf(mseAIC)
% ecdf(mseML)
% plot(range,Bcdf)
% plot(range,Fcdf)
% plot(range,AICcdf)
% xlabel('MSE')
% legend('Back','Forward','AIC','ML')

%% Statistical Properties
sympref('FloatingPointOutput',1);
format short
% Means
meanB = mean(mseB);
meanF = mean(mseF);
meanAIC = mean(mseAIC);
meanML = mean(mseML);
meanC = mean(mseC);
muVec = round([meanB; meanF; meanAIC; meanML; meanC],3);

% Variance
varB = var(mseB);
varF = var(mseF);
varAIC = var(mseAIC); 
varML = var(mseML);
varC = var(mseC);
varVec = round([varB; varF; varAIC; varML; varC],3);

% Medians
medB = median(mseB);
medF = median(mseF);
medAIC = median(mseAIC);
medML = median(mseML);
medC = median(mseC);
medVec = round([medB; medF; medAIC; medML; medC],3);



orderVec = round([mean(psaveF); mean(psaveF); mean(psaveAIC); mean(psaveML); mean(psaveC)],2);

TVarVec = round([mean(tvarB); mean(tvarF); mean(tvarAIC); mean(tvarML); mean(tvarC)],3);





comVec = [muVec, varVec, medVec, orderVec];
Criterion = {'BVIC (backcast)'; 'BVIC (forecast)'; 'AIC' ;'Maximum Likelihood'; 'AICc'};
T = table(Criterion, muVec, varVec, medVec, orderVec,TVarVec,'VariableNames',{'Criterion','Mean (MSE) ','Variance (MSE)' ,'Median (MSE)','Average order','Theoretical Variance'})
writetable(T,'setting1.csv')



