%% AUTHOR: Max Sibeijn. DATE: 03/05/21

clear all 
close all
clc


addpath('StructFiles')
addpath('Seizure_DATA')
load('clinical_metadata.mat')

rng default
tic
savestruct = 0; % 1 to save struct file, 0 to not save
%% SETTINGS

mm = 100; % number of steps used for backcasting
fm = 100; % number of steps used for forecasting
h_r = 1:fm;
h_range = h_r;

Windows = 20; %Number of windows
N = 1000;%nscale*(order+mm) + 2*mm; % change depending on size of windows
pmax = 80;%round(1.33*order);




Alpha = 100/(N-2*mm);
ratio = 2;
% BVIC parameters
Beta1 = 1*ceil(ratio*Alpha);
Beta2 = 100*ceil(ratio*Alpha);
Gamma = Alpha;

qrange = [1];
for q = qrange

load("MAYO010-interictal-block-"+q+".mat"); % load the specific dataset to be analysed
Ndata = Windows*N;%5000; %floor(length(evData)/1000)*1000; % Usable length of Ecog timeseries
Tstart = 10000;

MC = [5:8 13:16];%[17 18 19 20];%size(evData,1);
MonteCarloLoops = MC;
Windows = Ndata/N; % number of windows
    


Yu = evData(MC,1+Tstart:Ndata+Tstart);
Yu = Yu - mean(Yu,2);
%Notch Filter
wo = 60/(500/2);
bw = wo/35;
Ab = 3;
[fnum,fden] = iirnotch(wo,bw,Ab);
Y = filter(fnum,fden,Yu')';


%Normalize data
Y = (Y-mean(Y,2));%./std(Y,0,2);

  
  %create matrix for h
hsize = numel(h_range);
h_save1 = zeros(hsize,Windows*length(MC));
h_save2 = zeros(hsize,Windows*length(MC));
h_saveAIC = zeros(hsize,Windows*length(MC));
h_saveBIC = zeros(hsize,Windows*length(MC));
h_saveC = zeros(hsize,Windows*length(MC));


%% loop
counter = 0;
for mc = 1:length(MC)
   
    X = Y(mc,:);
   
    for ii = 1:Ndata/N
    XM(ii,:) = X((ii-1)*N+1:ii*N); 
    end

%% Model
for zz = 1:Windows

    Xloop = XM(zz,:)'; % The samples belonging to the zz'th window
     Xloop = (Xloop - mean(Xloop));%/std(Xloop);
    NT = N-fm;%ceil(0.5*(1.2*order+2*mm)+.5); % Determines the size of forecast set 
    NT2 = N-fm-mm;
    XT = zeros(2*NT-N,1);
    XT = Xloop(1:NT);
    XT2 = XT(mm+1:end);
    
   [pp_h(mc,zz), pp_p(mc,zz)] = pptest(Xloop);
   [lmc_h(mc,zz), lmc_p(mc,zz)] = lmctest(Xloop);
   [kpss_h(mc,zz), kpss_p(mc,zz)] = kpsstest(Xloop);
   
p = 1:pmax;%order+ceil(0.2*order);
idx = 1:numel(p);
acf = autocorr(XT,'NumLags',length(XT)-1);
acf2 = autocorr(XT2,'NumLags',length(XT2)-1);
% logL = zeros(length(p),1);
biasX = zeros(length(p),mm);
varX = zeros(length(p),mm);
        

for i = 1:numel(p)


    % LogLikelihood computation by OLS: 
    %Choose Xloop for in-sample forecast/ Choose XT for out-of-sample
    %forecast
    Xk = XT; 
    Xk2 = XT2;
    rL = acf(2:p(i)+1);
    XLS = hankel(flip(Xk(p(i):end-1)),flip(Xk(1:p(i)))); % tranposed hankel matrix
    YLS = flip(Xk(p(i)+1:end));
    phiLS = (XLS'*XLS)\XLS'*YLS;
     
     
    XLS2 = hankel(flip(Xk2(p(i):end-1)),flip(Xk2(1:p(i)))); % tranposed hankel matrix
    YLS2 = flip(Xk2(p(i)+1:end));
    phiLS2 = (XLS2'*XLS2)\XLS2'*YLS2;
%     phi_burg = arburg(Xk,i);
%     phiLS = -phi_burg(2:end)';
%    
%      phiLS = lsqr(XLS,YLS,1e-6,25);



    ISL = length(Xk)-p(i);             % in-sample size
    ISL2 = length(Xk2)-p(i);
    
%     eps = zeros(ISL,1);
    for k = 1:ISL
    eps(k,i) = Xk(p(i)+k) - phiLS'*flip(Xk(k:p(i)+k-1));
    end
    
    for k2 = 1:ISL2
    eps2(k2,i) = Xk2(p(i)+k2) - phiLS2'*flip(Xk2(k2:p(i)+k2-1));
    end
    
    sigma2(i) = (ISL)\sum(eps(:,i).^2);
    sigma2_2(i) = ISL2\sum(eps2(:,i).^2);

    logL(i,zz) = -length(XT)*log(sigma2(i)');
    logL2(i,zz) = -length(XT2)*log(sigma2_2(i)');


       R_n = toeplitz(acf2(1:p(i))); %AUTOCORRELATION MATRIX R
    for m = 1:mm
        rho_n = acf2(m+1:m+p(i));        %AUTOCORRELATION FUNCTION rho
        phi_f = R_n\rho_n;              % PHI autoregressive coefficient vector
        errorX(i,m) = Xloop(-m+mm+1) - phi_f'*XT2(1:p(i));
        varX(i,m)   = var(XT2)*(1- phi_f'*rho_n);
    end
    
    
    
    b2X(i,1) = norm(errorX(i,:))^2/mm;
    vX(i,1)  = sum(varX(i,:))/mm;
end


mlogL = logL(:,zz);
mlogL2 = logL2(:,zz);
mmse = mean(b2X,2);
mvar = mean(vX,2);
   
% end
    
% [mlogL, mmse, mvar, acf] = bcast(Xloop,XT,N,NT,mm,order);
    



%% Criterion

[~,pml] = max(mlogL);
[~,pml2] = max(mlogL2);
% [~,ptest] = min(mlogL);
% 
[~,pmse] = min(mmse);
[~,pvar] = min(mvar);




if mlogL(pml) > 0
    c = -1;
else
    c = 1;
end

% penal2 = mmse + (mmse.^2 + mmse);

fvec = [-mlogL2/abs(mlogL2(pml)), mmse/mmse(pml), mvar/mvar(pml)];
% fvec = [c*mlogL/mlogL(pml), penal2/penal2(pml), mvar/mvar(pml)];


W1 = [1; Beta1; Gamma];
W2 = [1; Beta2; Gamma];
[obj1, pmin1] = min(fvec*W1); %+p'/30 % BVIC
[obj2, pmin2] = min(fvec*W2); %+p'/30 % BVIC


aic2 = -mlogL + 2*p';
[aic, bic] = aicbic(mlogL/2,p',length(XT));
aicc = aic + (2*p'.^2+2*p')./(ISL - p' - 1);

[~,paic] = min(aic);
[~,pbic] = min(bic);
paic = p(paic);

[~,pc] = min(aicc);
pc = p(pc);




%% Forecast

pf = p(pmin1);
pf2 = p(pmin2);
pml = pbic;%p(pml);

%% YULE-WALKER METHOD
% [f_err, f_aic, f_ml, f_c, Va, Vaic, Vml, Vc] = fcast(Xloop,XT,acf,pf,paic,pml,pc,mm,NT);

for fh = 1:fm
    
    rho_A = acf2(fh+1:pf+fh);
    R_A = toeplitz(acf2(1:pf));
    phi_A = R_A\rho_A;
    
    rho_A2 = acf2(fh+1:pf2+fh);
    R_A2 = toeplitz(acf2(1:pf2));
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
    Va2(fh) = var(XT2)*(1 - phi_A2'*rho_A2);

    Xaic(fh) = phi_aic'*flip(Xloop(NT-paic+1:NT)); 
    Vaic(fh) = var(XT)*(1 - phi_aic'*rho_aic);
    
    Xml(fh) = phi_ml'*flip(Xloop(NT-pml+1:NT)); 
    Vml(fh) = var(XT)*(1 - phi_ml'*rho_ml);
    
    
    Xc(fh) = phi_c'*flip(Xloop(NT-pc+1:NT)); 
    Vc(fh) = var(XT)*(1 - phi_c'*rho_c);
    
end



Xf = Xloop(NT+h_range(1):NT+h_range(end));
f_err = Xa(h_range)'- Xf;
f_err2 = Xa2(h_range)'- Xf;
f_aic = Xaic(h_range)'-Xf;
f_ml = Xml(h_range)'-Xf;
f_c = Xc(h_range)'-Xf;

% effective errors
Norm = mean(abs(diff(XT2)));
f_err = f_err/Norm;
f_err2 = f_err2/Norm;
f_aic = f_aic/Norm;
f_ml = f_ml/Norm;
f_c = f_c/Norm;

h_save1(:,zz+Windows*(mc-1))= abs(f_err);
h_save2(:,zz+Windows*(mc-1))= abs(f_err2);
h_saveAIC(:,zz+Windows*(mc-1))= abs(f_aic);
h_saveBIC(:,zz+Windows*(mc-1))= abs(f_ml);
h_saveC(:,zz+Windows*(mc-1))= abs(f_c);




%% Validation test
MSE_F(zz,mc) = mean(f_err.^2);
MSE_F2(zz,mc) = mean(f_err2.^2);

MSE_B(zz,mc) = mmse(pmin1);
MSE_AIC(zz,mc) = mean(f_aic.^2);
MSE_ML(zz,mc) = mean(f_ml.^2);
MSE_C(zz,mc) = mean(f_c.^2);

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

mh_1(:,mc) = mean(h_save1(:,Windows*(mc-1)+1:Windows*mc),2);
mh_2(:,mc) = mean(h_save2(:,Windows*(mc-1)+1:Windows*mc),2);
mh_3(:,mc) = mean(h_saveAIC(:,Windows*(mc-1)+1:Windows*mc),2);
mh_4(:,mc) = mean(h_saveBIC(:,Windows*(mc-1)+1:Windows*mc),2);
mh_5(:,mc) = mean(h_saveC(:,Windows*(mc-1)+1:Windows*mc),2);

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

[pt(mc),h1(mc)] = signrank(MSE_F(:,mc),MSE_AIC(:,mc));
[pt2(mc),h2(mc)] = signrank(MSE_F(:,mc),MSE_ML(:,mc));
[pt3(mc),h3(mc)] = signrank(MSE_F(:,mc),MSE_C(:,mc));
[pt4(mc),h4(mc)] = signrank(MSE_F(:,mc),MSE_F2(:,mc));

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

hVec = [NaN; sum(h4); sum(h1); sum(h2); sum(h3)]*100;



comVec = [muVec, varVec, TVarVec, orderVec];
Criterion = {"BVIC (beta = "+Beta1+")"; "BVIC (beta ="+Beta2+")"; 'AIC' ;'BIC'; 'AICc'};
T = table(Criterion, muVec, varVec, medVec, orderVec,TVarVec, hVec, 'VariableNames',{'Criterion','Mean (MSE) ','Variance (MSE)' ,'Median (MSE)','Av. order','Th. Variance','h-Values'});
% writetable(T,'setting1.csv')

% 
% for k=1:4
%     fprintf('%8.3f & %8.2f & %8.2f & %8.1f  \n',comVec(k+1,:))
% end


% figure
% plot(Xloop)
% hold on
% plot(N-NT+1:NT,XT,'LineWidth',1.5)
% plot(NT+1:NT+fm,Xa,'LineWidth',1.5)
% plot(NT+1:NT+fm,Xc,'LineWidth',1.5)
% legend('X_{window}','X_{fit}','F_{bvic}','F_{aicc}')


%% create struct

if savestruct > 0.5


INFO = struct('SaveStructure','Info-ECoG-order-horizon-pmax-N.mat','Horizon',mm,'samplesize',N,'Beta1',Beta1,'Beta2',Beta2,'Gamma',Gamma,'MSE_BVIC',MSE_F,'MSE_AIC',MSE_AIC,'MSE_ML',MSE_ML,'MSE_C',MSE_C,'pmax',pmax,'TABLE',T);
filename = "StructFiles/h_ranges/ECoG-start-" + h_range(1) + "-h-"+ mm + "-pm-" + pmax + "-N-" + N +".mat";
save(filename,'INFO');
end




fvec(:,1) = fvec(:,1) + 2;

h_bvic(:,length(MC)*(q-1)+1:length(MC)*q) = mh_1;
h_aic(:,length(MC)*(q-1)+1:length(MC)*q) = mh_3;
end

%% PLOT
MSFE = [mean(h_save1,2), mean(h_save2,2), mean(h_saveAIC,2), mean(h_saveBIC,2), mean(h_saveC,2)];

Mb = mean(mh_1,2);
Ma = mean(mh_3,2);
Varb = std(mh_1')';
Vara = std(mh_3')';
% % 
addpath('ECoGplots')
%save("ECoGplots/mase_MAYO016-interictal-"+qrange(1)+"-filtered.mat",'mh_1','mh_2','mh_3','mh_4','mh_5','Y','Xloop','Xa','Va','Tstart')
% save("ECoGplots/mase_HUP86-ictal-"+qrange(1)+"onechannel.mat",'h_save1')
figure


plot(h_r,Mb,'b','LineWidth',2)
hold on

% plot(h_r,Ma,'r','LineWidth',1)
% plot(h_r,Mb+Varb,'--b',h_r,Mb-Varb,'--b')

% plot(h_r,Ma,'r','LineWidth',2)
% hold on
% plot(h_r,Ma+1.96*Vara,'--r',h_r,Ma-1.96*Vara,'--r')



hold on

CBf_top = Mb+1.96*Varb;
CBf_bot = Mb-1.96*Varb;

x2 = [h_r, flip(h_r)];
inBetween = [CBf_top', flip(CBf_bot')];

pat1 = patch(x2,inBetween,[.85 0.85 1],'EdgeColor',[1 1 1])
uistack(pat1,'bottom')




 legend('95$\%$ confidence (1.96$\sigma$)','Average MSFE for BVIC','Interpreter','latex','FontSize',12)
plots=get(gca, 'Children');
legend(plots([1, 2]), {'Sample mean (BVIC)', 'Sample deviation (95$\%$ confidence)'});


% ylim([0 0.35])
xlabel('Forecasted step (h)')
ylabel('MASE (-)')
% title("MAYO016-BLOCK-"+qrange(1)+" (ictal)")

MSFEn = MSFE./MSFE(:,3);

% figure
% plot(h_r, MSFE','LineWidth',1.5)
% xlabel('h')
% ylabel('Average MSFE')
% title('Average MSFE for each prediction horizon')
% legend("BVIC(\beta="+Beta1+")","BVIC(\beta="+Beta2+")",'AIC','BIC','AICc')
% 
% 
% 
% figure
% plot(h_r, MSFEn','LineWidth',1.5)
% xlabel('h')
% ylabel('Average MSFE')
% title('Average MSFE for each prediction horizon (normalized)')
% legend("BVIC(\beta="+Beta1+")","BVIC(\beta="+Beta2+")",'AIC','BIC','AICc')

%% 

figure
plot(Xloop)
hold on
plot(NT+1:NT+fm,Xa);
plot(NT+1:NT+fm,Xa+CB,'--r',NT+1:NT+fm,Xa-CB,'--r')

% 
% wo = 60/(500/2);
% bw = wo/20;
% [fnum,fden] = iirnotch(wo,bw);
% Yfilt = filter(fnum,fden,Y(end,:));
% figure
% pspectrum(Yfilt)
% kpss = sum(kpss_h(:))/length(kpss_h(:))
% 
% lmc = sum(lmc_h(:))/length(lmc_h(:))
% ppt = sum(pp_h(:))/length(pp_h(:))