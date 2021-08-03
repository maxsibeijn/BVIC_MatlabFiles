clear all
close all
clc

addpath('Seizure_DATA')
% 
rng default

DataType = 1; % 0 if generated(small order), 1 if generated(big order), 2 if real data

%% SETTINGS
order = 25; %Order of the generated model
mm = 5; % number of steps used for forecasting
nscale = 2; %size of training set = nscale*(p+h)
alpha = 0; % noise parameter

% BVIC parameters
Beta1 = 3;
Beta2 = 50;
Gamma = 1;


%% ALGORITHM

if DataType == 0
    N = 100; % size of window
    MonteCarloLoops = 300;
    Windows = 100; % number of windows
elseif DataType == 1
    N = nscale*(order+mm)+2*mm;%200; % size of window
    MonteCarloLoops = 100;
    Windows = 1; % number of windows
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
[X,XM,phi,pol,A] = TSgen(N,order, Windows, alpha);

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

    Xloop = XM(:,zz); % The 100 samples belonging to the zz'th window
    Xloop = (Xloop - mean(Xloop))/std(Xloop);
    NT = 0.5*N+round(nscale/2*(order+mm));%ceil(0.5*(1.2*order+2*mm)+.5); % Determines the size of forecast set 
    XT = zeros(2*NT-N,1);
    XT = Xloop(N-NT+1:NT);
    XT2 = XT(mm+1:end);
    
    
    
% function [mlogL, mmse, mvar, acf] = bcast(Xloop,XT,N,NT,mm,order)


p = 1:order+ceil(0.2*order);
idx = 1:numel(p);
acf = autocorr(XT,'NumLags',max(p)+mm);
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
    Xk = XT; 
    rL = acf(2:p(i)+1);

% %      XLSc = toeplitz(Xk(p(i):end-1),flip(Xk(1:p(i))));
     XLS = hankel(flip(Xk(p(i):end-1)),flip(Xk(1:p(i)))); % tranposed hankel matrix
     YLS = flip(Xk(p(i)+1:end));
     phiLS = (XLS'*XLS)\XLS'*YLS;
%     phi_L = R_n\rL;
%     phiLS = arburg(XT,p(i));
%     phiLS = phiLS(2:end)';
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
       R_n = toeplitz(acf2(1:p(i))); %AUTOCORRELATION MATRIX R
    for m = 1:mm
        rho_n = acf2(m+1:m+p(i));        %AUTOCORRELATION FUNCTION rho
        phi_f = R_n\rho_n;              % PHI autoregressive coefficient vector
        errorX(i,m) = Xloop(N-NT-m+mm+1) - phi_f'*XT2(1:p(i));
        varX(i,m)   = var(XT2)*(1- phi_f'*rho_n);
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


if mlogL(pml) > 0
    c = -1;
else
    c = 1;
end

fvec = [-mlogL/abs(mlogL(pml)), mmse/mmse(pml), mvar/mvar(pml)];




W1 = [1; Beta1; Gamma];
W2 = [1; Beta2; Gamma];
[obj1, pmin1] = min(fvec*W1); %+p'/30 % BVIC
[obj2, pmin2] = min(fvec*W2); %+p'/30 % BVIC

aic2 = -mlogL + 2*p';
[aic, bic] = aicbic(mlogL/2,p',ISL);
aicc = aic + (2*p'.^2+2*p')./(ISL - p' - 1);

[~,paic] = min(aic);
[~,pbic] = min(bic);
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

pf = p(pmin1);
pf2 = p(pmin2);

pml = pbic;%p(pml);

%% YULE-WALKER METHOD
% [f_err, f_aic, f_ml, f_c, Va, Vaic, Vml, Vc] = fcast(Xloop,XT,acf,pf,paic,pml,pc,mm,NT);
fm = mm;
for fh = 1:fm
    
    rho_A = acf(fh+1:pf+fh);
    R_A = toeplitz(acf(1:pf));
    phi_A = R_A\rho_A;
    
        
    rho_A2 = acf(fh+1:pf2+fh);
    R_A2 = toeplitz(acf(1:pf2));
    phi_A2 = R_A2\rho_A2;
    
    rho_aic = acf(fh+1:paic+fh);
    R_aic = toeplitz(acf(1:paic));
    phi_aic = R_aic\rho_aic;
    
    
    rho_ml = acf(fh+1:pml+fh);
    R_ml = toeplitz(acf(1:pml));
    phi_ml = R_ml\rho_ml;
    
    rho_c = acf(fh+1:pc+fh);
    R_c = toeplitz(acf(1:pc));
    phi_c = R_c\rho_c;
    
    Xa(fh) = phi_A'*flip(Xloop(NT-pf+1:NT)); 
    Va(fh) = var(XT)*(1 - phi_A'*rho_A);
    
    Xa2(fh) = phi_A2'*flip(Xloop(NT-pf2+1:NT)); 
    Va2(fh) = var(XT)*(1 - phi_A2'*rho_A2);


    Xaic(fh) = phi_aic'*flip(Xloop(NT-paic+1:NT)); 
    Vaic(fh) = var(XT)*(1 - phi_aic'*rho_aic);
    
    Xml(fh) = phi_ml'*flip(Xloop(NT-pml+1:NT)); 
    Vml(fh) = var(XT)*(1 - phi_ml'*rho_ml);
    
    
    Xc(fh) = phi_c'*flip(Xloop(NT-pc+1:NT)); 
    Vc(fh) = var(XT)*(1 - phi_c'*rho_c);
    
end


Xf = Xloop(NT+1:NT+fm);
f_err = Xa'- Xf;
f_err2 = Xa2'- Xf;
f_aic = Xaic'-Xf;
f_ml = Xml'-Xf;
f_c = Xc'-Xf;
%% BURG METHOD
% fm = mm;
% yule = zeros(fm,paic+fm+1);
% levi = zeros(fm,paic+fm+1);
% burg = zeros(fm,paic+fm+1);
% araic = zeros(fm,paic);
% 
% burgf = arburg(XT,pf);
% burgf = -burgf(2:end);
% 
% [burgaic, ee] = arburg(XT,paic);
% burgaic = -burgaic(2:end);
% 
% burgc = arburg(XT,pc);
% burgc = -burgc(2:end);
% 
% burgml = arburg(XT,pml);
% burgml = -burgml(2:end);
% 
% 
% burgf_save = burgf;
% burgf_update = burgf;
% 
% burgaic_save = burgaic;
% burgaic_update = burgaic;
% 
% burgc_save = burgc;
% burgc_update = burgc;
% 
% burgml_save = burgml;
% burgml_update = burgml;
% yulephi_save = yule;
% yulephi_update = yule;
% 
% 
% for fh = 1:fm
%     
%     rho_A = acf(fh+1:pf+fh);
%     R_A = toeplitz(acf(1:pf));
%     phi_A = R_A\rho_A;
%     
%     rho_aic = acf(fh+1:paic+fh);
%     R_aic = toeplitz(acf(1:paic));
%     phi_aic = R_aic\rho_aic;
%     
%     araic(fh,:) = phi_aic';
% 
%     rho_ml = acf(fh+1:pml+fh);
%     R_ml = toeplitz(acf(1:pml));
%     phi_ml = R_ml\rho_ml;
%     
%     rho_c = acf(fh+1:pc+fh);
%     R_c = toeplitz(acf(1:pc));
%     phi_c = R_c\rho_c;
%     
%     Xa(fh) = burgf_update*flip(Xloop(NT-pf+1:NT)); 
%     Va(fh) = var(XT)*(1 - burgf_update*rho_A);
% 
%     Xaic(fh) = burgaic_update*flip(Xloop(NT-paic+1:NT)); 
%     Vaic(fh) = var(XT)*(1 - burgaic_update*rho_aic);
%     
%     Xml(fh) = burgml_update*flip(Xloop(NT-pml+1:NT)); 
%     Vml(fh) = var(XT)*(1 - burgml_update*rho_ml);
%     
%     
%     Xc(fh) = burgc_update*flip(Xloop(NT-pc+1:NT)); 
%     Vc(fh) = var(XT)*(1 - burgc_update*rho_c);
%     
%     burgf_update = burgf*burgf_update(1) + [burgf_update(2:end) 0];
%     burgf_save = [burgf_save; burgf_update];
%     
%     burgaic_update = burgaic*burgaic_update(1) + [burgaic_update(2:end) 0];
%     burgaic_save = [burgaic_save; burgaic_update];
%     
%     burgc_update = burgc*burgc_update(1) + [burgc_update(2:end) 0];
%     burgc_save = [burgc_save; burgc_update];
%     
%     burgml_update = burgml*burgml_update(1) + [burgml_update(2:end) 0];
%     burgml_save = [burgml_save; burgml_update];
%     
% 
% end
% 
% Xf = Xloop(NT+1:NT+fm);
% f_err = Xa'- Xf;
% f_aic = Xaic'-Xf;
% f_ml = Xml'-Xf;
% f_c = Xc'-Xf;

%% Validation test
MSE_F(zz,mc) = norm(f_err)^2/mm;
MSE_F2(zz,mc) = norm(f_err2)^2/mm;

MSE_B(zz,mc) = mmse(pmin1);
MSE_AIC(zz,mc) = norm(f_aic)^2/mm;
MSE_ML(zz,mc) = norm(f_ml)^2/mm;
MSE_C(zz,mc) = norm(f_c)^2/mm;

V_F(zz,mc) = mean(Va); 
V_F2(zz,mc) = mean(Va2); 
V_B(zz,mc) = mvar(pmin1);
V_AIC(zz,mc) = mean(Vaic);
V_ML(zz,mc) = mean(Vml);
V_C(zz,mc) = mean(Vc);

CB = 1.96*sqrt(Va);


p1(zz,mc) = pf;
p12(zz,mc) = pf2;
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
MaxC(mc) = max(MSE_C(:,mc));

MinB(mc) = min(MSE_B(:,mc));
MinF(mc) = min(MSE_F(:,mc));
MinAIC(mc) = min(MSE_AIC(:,mc));
MinML(mc) = min(MSE_ML(:,mc));
MinC(mc) = min(MSE_C(:,mc));


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


[pt1(mc),h1] = signrank(MSE_F,MSE_AIC);
[pt2(mc),h2] = signrank(MSE_F,MSE_ML);
[pt3(mc),h3] = signrank(MSE_F,MSE_C);
[pt4(mc),h4] = signrank(MSE_F,MSE_F2);

% One vector containng all MSE estimates
mseB = MSE_B(:);
mseF = MSE_F(:);
mseF2 = MSE_F2(:);
mseAIC = MSE_AIC(:);
mseML = MSE_ML(:);
mseC = MSE_C(:);

% vector containing all theoretical variance estimates
tvarB = V_B(:);
tvarF = V_F(:);
tvarF2 = V_F2(:);
tvarAIC = V_AIC(:);
tvarML = V_ML(:);
tvarC = V_C(:);



psaveF = p1(:);
psaveF2 = p12(:);
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
% % 
% % 
% figure
% BW2 = .1;
% M = max(MaxAIC);
% ax = [0 M+.5 0 15];
% hline = median(MaxF);
% haic = median(MaxAIC);
% hbic = median(MaxML);
% hC = median(MaxC);
% hold on
% % subplot(4,1,1), histogram(MaxB,'BinWidth',BW2)
% % axis(ax)
% 
% title('BVIC (backcast)')
% 
% subplot(4,1,1), 
% xline(hline, '--k','LineWidth',2)
% hold on
% histogram(MaxF,'BinWidth',BW2)
%   axis(ax)
%   legend('Median')
% title('BVIC')
% 
% subplot(4,1,2),
% xline(haic, '--k','LineWidth',2)
% hold on
%  histogram(MaxAIC,'BinWidth',BW2)
% axis(ax)
% 
% legend('Median')
% title('AIC')
% 
% subplot(4,1,3),
% xline(hbic, '--k','LineWidth',2)
% hold on
%  histogram(MaxML,'BinWidth',BW2)
% axis(ax)
% legend('Median')
% title('BIC')
% 
% subplot(4,1,4), xline(hC, '--k','LineWidth',2)
% hold on
% histogram(MaxC,'BinWidth',BW2)
% 
% axis(ax)
% xlabel('MSE')
% title('AICc')
% legend('Median')
% % % 
% % % 
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
meanF2 = mean(mseF2);
meanAIC = mean(mseAIC);
meanML = mean(mseML);
meanC = mean(mseC);
muVec = round([meanF; meanF2; meanAIC; meanML; meanC],3);

% Variance
varB = var(mseB);
varF = var(mseF);
varF2 = var(mseF2);
varAIC = var(mseAIC); 
varML = var(mseML);
varC = var(mseC);
varVec = round([varF; varF2; varAIC; varML; varC],3);

% Medians
medB = median(mseB);
medF = median(mseF);
medF2 = median(mseF2);
medAIC = median(mseAIC);
medML = median(mseML);
medC = median(mseC);
medVec = round([medF; medF2; medAIC; medML; medC],3);



orderVec = round([mean(psaveF); mean(psaveF2); mean(psaveAIC); mean(psaveML); mean(psaveC)],2);

TVarVec = round([mean(tvarF); mean(tvarF2); mean(tvarAIC); mean(tvarML); mean(tvarC)],3);

hVec = [NaN; sum(h4); sum(h1); sum(h2); sum(h3)]*100/MonteCarloLoops;


comVec = [muVec, varVec, TVarVec, orderVec];
Criterion = {'BVIC (beta = 1)'; 'BVIC (beta = 5)'; 'AIC' ;'BIC'; 'AICc'};
T = table(Criterion, muVec, varVec, medVec, orderVec,TVarVec, hVec, 'VariableNames',{'Criterion','Mean (MSE) ','Variance (MSE)' ,'Median (MSE)','Av. order','Th. Variance','h-Values'})
% writetable(T,'setting1.csv')



fprintf('%8.3f & %8.2f & %8.2f & %8.1f  \n',comVec')



figure
plot(Xloop)
hold on
plot(N-NT+1:NT,XT,'LineWidth',1.5)
plot(NT+1:NT+fm,Xa,'LineWidth',1.5)
plot(NT+1:NT+fm,Xc,'LineWidth',1.5)
legend('X_{window}','X_{fit}','F_{bvic}','F_{aicc}')


%t-test
% [h,pt] = ttest2(mseF,mseAIC);

% delta_dB = -round(10*log10(alpha));
% INFO = struct('SaveStructure','Info-exp-poles-delta-beta-gamma.mat','Order',order,'Horizon',mm,'delta',alpha,'Training_set',N,'Poles',pol,'Beta',Beta,'Gamma',Gamma,'MSE_BVIC',MSE_F,'MSE_AIC',MSE_AIC,'MSE_ML',MSE_ML,'MSE_C',MSE_C,'TABLE',T);
% filename = "StructFiles/Info-E1-P4-" +  delta_dB + "-" + Beta + "-" + Gamma + ".mat";
% save(filename,'INFO');
        