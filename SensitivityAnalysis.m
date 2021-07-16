%% AUTHOR: Max Sibeijn. DATE: 31/03/21

clear all 
close all
clc

addpath('StructFiles')

rng default
tic
savestruct = 0; % 1 to save struct file, 0 to not save
%% SETTINGS

rng default
order = 30; %Order of the generated model
mm = ceil(1*order); % number of steps used for forecasting
h_r = 1:mm;
h_range = h_r;
nscale = 4.5; %size of training set = nscale*(p+h)
delta = 0.1; % noise parameter
Windows = 100; %Number of windows
N = ceil(nscale*(order*2))+2*mm; % change depending on size of windows
pmax = round(1.2*order);
MC = 100;
MonteCarloLoops = MC;

Alpha = 100/(N-2*mm);
ratio = 2;
% BVIC parameters
Beta1 = 1;%*round(ratio*Alpha);
Beta2 = 5;%*round(ratio*Alpha);
Gamma = 1;%Alpha;

P = zeros(Windows,order);
Poly = zeros(Windows,order+1);

Nts = N*MC; % length of X

X = zeros(Windows,3*Nts);
Xnorm = zeros(Windows,3*Nts);
Y = zeros(Windows,3*Nts);

%create matrix for h
hsize = numel(h_range);
h_save1 = zeros(hsize,Windows*MC);
h_save2 = zeros(hsize,Windows*MC);
h_saveAIC = zeros(hsize,Windows*MC);
h_saveBIC = zeros(hsize,Windows*MC);
h_saveC = zeros(hsize,Windows*MC);

for j = 1:Windows
    for k = 1:order/2
        mag = 0.5+0.49*rand; % magnitude of the pole
        angle = pi*rand;%2*pi*(k-1)/order; % phase angle of the pole
        c_real = cos(angle)*mag;
        c_imag = sin(angle)*mag;
        
        if rand > 0.5
            P(j,2*k - 1) = c_real + c_imag*1i;
            P(j,2*k) = c_real - c_imag*1i;
        else
            P(j,2*k - 1) = sign(rand-0.5)*(rand*0.5+0.49);%c_real;
            P(j,2*k) = sign(rand-0.5)*(rand*0.5+0.49);%-c_real;
%             P(j,2*k - 1) = sign(rand-0.5)*(0.1+rand*0.75)+0i;
%             P(j,2*k) = sign(rand-0.5)*(0.1+rand*0.75)+0i;
        end
    end
    W = normrnd(0,1,[1 3*Nts]); % Generate error dynamics
    Z = normrnd(0,1,[1 3*Nts]); % measurement noise
    Poly(j,:) = poly(P(j,:));
    X(j,:) = filter(1,Poly(j,:),W);
    Xnorm(j,:) = (X(j,:)-mean(X(j,:)))/std(X(j,:));
    Y(j,:) = Xnorm(j,:) + delta*Z;
end


Y = Y(:,2*Nts+1:end);
[rows,~] = find(isnan(Y));
nanrows = unique(rows);
Y(nanrows,:) = [];%Y(nanrows-1,:);

plot(mean(abs(Poly)))

%%
counter = 0;
for mc = 1:MC
    XM = Y(:,1+N*(mc-1):N*mc);
    XM = XM';



%% Model
for zz = 1:Windows-length(nanrows)

    Xloop = XM(:,zz); % The samples belonging to the zz'th window
    Xloop = (Xloop - mean(Xloop));%/std(Xloop);
    NT = 0.5*N+(N/2-mm);%ceil(0.5*(1.2*order+2*mm)+.5); % Determines the size of forecast set 
    XT = zeros(2*NT-N,1);
    XT = Xloop(N-NT+1:NT);
    XT2 = XT(mm+1:end);
    
   


p = 1:pmax;%order+ceil(0.33*order);
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



%     ISL = length(Xk)-max(p);             % in-sample size
    ISL = length(Xk)-p(i);%max(p); 
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
        errorX(i,m) = Xloop(N-NT-m+mm+1) - phi_f'*XT2(1:p(i));
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





% penal2 = mmse + (mmse.^2 + mmse);

fvec = [-mlogL2/abs(mlogL2(pml2)), mmse/mmse(pml2), mvar/mvar(pml2)];
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


%% YULE-WALKER METHOD
% [f_err, f_aic, f_ml, f_c, Va, Vaic, Vml, Vc] = fcast(Xloop,XT,acf,pf,paic,pml,pc,mm,NT);
fm = mm;
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
    
    
    rho_ml = acf(fh+1:pbic+fh);
    R_ml = toeplitz(acf(1:pbic));
    phi_ml = R_ml\rho_ml;
    
    rho_c = acf(fh+1:pc+fh);
    R_c = toeplitz(acf(1:pc));
    phi_c = R_c\rho_c;
    
    Xa(fh) = phi_A'*flip(Xloop(NT-pf+1:NT)); 
    Va(fh) = var(XT2)*(1 - phi_A'*rho_A);
    
    Xa2(fh) = phi_A2'*flip(Xloop(NT-pf2+1:NT)); 
    Va2(fh) = var(XT2)*(1 - phi_A2'*rho_A2);

    Xaic(fh) = phi_aic'*flip(Xloop(NT-paic+1:NT)); 
    Vaic(fh) = var(XT)*(1 - phi_aic'*rho_aic);
    
    Xml(fh) = phi_ml'*flip(Xloop(NT-pbic+1:NT)); 
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

Nrm = mean(abs(diff(XT2)));
h_save1(:,zz+Windows*(mc-1))= (f_err).^2;
h_save2(:,zz+Windows*(mc-1))= (f_err2).^2;
h_saveAIC(:,zz+Windows*(mc-1))= (f_aic).^2;
h_saveBIC(:,zz+Windows*(mc-1))= (f_ml).^2;
h_saveC(:,zz+Windows*(mc-1))= (f_c).^2;

% 
% h_save1(:,mc+MC*(zz-1))= abs(f_err)/Nrm;
% h_save2(:,mc+MC*(zz-1))= abs(f_err2)/Nrm;
% h_saveAIC(:,mc+MC*(zz-1))= abs(f_aic)/Nrm;
% h_saveBIC(:,mc+MC*(zz-1))= abs(f_ml)/Nrm;
% h_saveC(:,mc+MC*(zz-1))= abs(f_c)/Nrm;

%% Validation test
MSE_F(zz,mc) = mean(f_err.^2);
MSE_F2(zz,mc) = mean(f_err2.^2);

MSE_B(zz,mc) = mmse(pmin1);
MSE_AIC(zz,mc) = mean(f_aic.^2);
MSE_ML(zz,mc) = mean(f_ml.^2);
MSE_C(zz,mc) = mean(f_c.^2);
% MSE_F(zz,mc) = mean(abs(f_err)/Nrm);
% MSE_F2(zz,mc) = mean(abs(f_err2)/Nrm);
% 
% MSE_B(zz,mc) = mmse(pmin1);
% MSE_AIC(zz,mc) = mean(abs(f_aic)/Nrm);
% MSE_ML(zz,mc) = mean(abs(f_ml)/Nrm);
% MSE_C(zz,mc) = mean(abs(f_c)/Nrm);


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
p3(zz,mc) = pbic;
p4(zz,mc) = pc;
end



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

gr = {'bvic1','bvic2','aic','bic','aicc'};
pKW(mc) = kruskalwallis([MSE_F(:,mc) MSE_F2(:,mc) MSE_AIC(:,mc) MSE_ML(:,mc) MSE_C(:,mc)],gr,'off');
pAN(mc) = anova1([MSE_F(:,mc) MSE_F2(:,mc) MSE_AIC(:,mc) MSE_ML(:,mc) MSE_C(:,mc)],gr,'off');
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

hVec = [NaN; sum(h4); sum(h1); sum(h2); sum(h3)]*100/MC;



comVec = [muVec, varVec, TVarVec, orderVec];
Criterion = {'BVIC (beta1)'; 'BVIC (beta2)'; 'AIC' ;'BIC'; 'AICc'};
T = table(Criterion, muVec, varVec, medVec, orderVec,TVarVec, hVec, 'VariableNames',{'Criterion','Mean (MSE) ','Variance (MSE)' ,'Median (MSE)','Av. order','Th. Variance','h-Values'});
% writetable(T,'setting1.csv')


fprintf("$BVIC (\\beta = "+Beta1+")$ & %8.3f & %8.2f & %8.2f & %8.1f \\\\ \n", comVec(1,:)')
fprintf("$BVIC (\\beta = "+Beta2+")$ & %8.3f & %8.2f & %8.2f & %8.1f \\\\  \n", comVec(2,:)')
fprintf('$AIC$ & %8.3f & %8.2f & %8.2f & %8.1f \\\\  \n', comVec(3,:)')
fprintf('$BIC$ & %8.3f & %8.2f & %8.2f & %8.1f \\\\ \n', comVec(4,:)')
fprintf('$AICc$ & %8.3f & %8.2f & %8.2f & %8.1f \\\\  \n', comVec(5,:)')




plot(Xloop)
hold on
plot(N-NT+1:NT,XT,'LineWidth',1.5)
plot(NT+1:NT+fm,Xa,'LineWidth',1.5)
plot(NT+1:NT+fm,Xc,'LineWidth',1.5)
legend('X_{window}','X_{fit}','F_{bvic}','F_{aicc}')

%% create struct

if savestruct > 0.5

delta_dB = -10*log10(delta);
INFO = struct('SaveStructure','Info-N-order-horizon-beta1-beta2-gamma.mat','Order',order,'Horizon',mm,'delta',delta,'Training_set',N,'Poles',P,'AR_Parameters',Poly,'Beta1',Beta1,'Beta2',Beta2,'Gamma',Gamma,'MSE_BVIC',MSE_F,'MSE_AIC',MSE_AIC,'MSE_ML',MSE_ML,'MSE_C',MSE_C,'p_select1',psaveF,'p_select2',psaveF2,'p_aicc',psaveC,'TABLE',T,'pANOVA',pAN,'pKW',pKW);
filename = "StructFiles/EXPERIMENT2/"+N+"-"+pmax+"-"+ order+ "-"+ mm +"-"+nscale+".mat";
save(filename,'INFO');
end

% if savestruct > 0.5
% 
% delta_dB = -10*log10(delta);
% INFO = struct('SaveStructure','Info-N-order-horizon-beta1-beta2-gamma.mat','Order',order,'Horizon',mm,'delta',delta,'Training_set',N,'Poles',P,'AR_Parameters',Poly,'Beta1',Beta1,'Beta2',Beta2,'Gamma',Gamma,'MSE_BVIC',MSE_F,'MSE_AIC',MSE_AIC,'MSE_ML',MSE_ML,'MSE_C',MSE_C,'p_select1',psaveF,'p_select2',psaveF2,'p_aicc',psaveC,'TABLE',T);
% filename =  "StructFiles/h_ranges/Synth-order-" + h_range(1) + "-h-"+ mm + "-pm-" + pmax + "-N-" + N +".mat";
% save(filename,'INFO');
% end

% MSFE(:,h_range) = muVec;
% 
% fvec(:,1) = fvec(:,1) + 2;
% 
% 
MSFE = [mean(h_save1,2), mean(h_save2,2), mean(h_saveAIC,2), mean(h_saveBIC,2), mean(h_saveC,2)] ;

MSFEn = MSFE./MSFE(:,3);

figure
plot(h_r, MSFE','LineWidth',1.5)
xlabel('h')
ylabel('Average MSFE')
title('Average MSFE for each prediction horizon')
legend("BVIC(r="+ratio+")","BVIC(r="+Beta2/Beta1*ratio+")",'AIC','BIC','AICc')


figure
plot(h_r, MSFEn','LineWidth',1.5)
xlabel('h')
ylabel('Average MSFE')
title('Average MSFE for each prediction horizon (normalized)')
legend("BVIC(r="+ratio+")","BVIC(r="+Beta2/Beta1*ratio+")",'AIC','BIC','AICc')

%% EXTRA
 % Histogram for maximum MSE

% 
% % 
% % 
% figure
% BW2 = .1;
% M = max([MaxAIC,MaxF,MaxML,MaxC]);
% ax = [0 M+.5 0 15];
% hline = mean(MaxF);
% haic = mean(MaxAIC);
% hbic = mean(MaxML);
% hC = mean(MaxC);
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
%   legend('Mean')
% title('BVIC')
% 
% subplot(4,1,2),
% xline(haic, '--k','LineWidth',2)
% hold on
%  histogram(MaxAIC,'BinWidth',BW2)
% axis(ax)
% 
% legend('Mean')
% title('AIC')
% 
% subplot(4,1,3),
% xline(hbic, '--k','LineWidth',2)
% hold on
%  histogram(MaxML,'BinWidth',BW2)
% axis(ax)
% legend('Mean')
% title('BIC')
% 
% subplot(4,1,4), xline(hC, '--k','LineWidth',2)
% hold on
% histogram(MaxC,'BinWidth',BW2)
% 
% axis(ax)
% xlabel('MSE')
% title('AICc')
% legend('Mean')




figure
histogram(psaveF)
title("\bf{BVIC($\beta$ = "+Beta1+", $\gamma$ = "+Gamma+")}",'FontSize',20,'Interpreter','latex')
xlabel('selected order','FontSize',16)
ylabel('number of times selected','FontSize',16)
 ylim([0 1200])
 set(gca,'FontSize',15);

figure
histogram(psaveF2)
title("\bf{BVIC($\beta$ = "+Beta2+", $\gamma$ = "+Gamma+")}",'FontSize',20,'Interpreter','latex')
xlabel('selected order','FontSize',16)
ylabel('number of times selected','FontSize',16)
ylim([0 1200])
set(gca,'FontSize',15);
% figure
% histogram(psaveML)
% title('BIC')
% xlabel('selected order by method')
figure
histogram(psaveC)
title("\bf{AICc}",'FontSize',20,'Interpreter','latex')
xlabel('selected order','FontSize',16)
ylabel('number of times selected','FontSize',16)
ylim([0 1200])
set(gca,'FontSize',15);
% mp1 = mean(p1,1);
% mp2 = mean(p2,1);
% figure
% histogram(mp1,20)
% xlim([0 15])
% title('BVIC')
% figure
% histogram(mp2,20)
% xlim([0 15])
% title('AICc')


% figure
% plot(fvec,'LineWidth',1.3)
% hold on
% plot(2*p/ISL,'LineWidth',1.3)
% plot((2*p'+ (2*p'.^2+2*p')./(ISL - p' - 1))/ISL,'LineWidth',1.3)
% xlabel('model order')
% ylabel('penalty magnitude')
% 
% legend('-logL','mse','var','AIC','AICc')
% title('Penalty terms')

%% kruskall wallis test
% gv = [mseF mseF2 mseAIC mseML mseC];
% ppkw = kruskalwallis(gv);
pkw = zeros(100,mm);
for ki = 1:mm
HH1 = h_save1(ki,:);
HH2 = h_save2(ki,:);
HH3 = h_saveAIC(ki,:);
HH4 = h_saveBIC(ki,:);
HH5 = h_saveC(ki,:);


H1 = reshape(HH1,100,100);
H2 = reshape(HH2,100,100);
H3 = reshape(HH3,100,100);
H4 = reshape(HH4,100,100);
H5 = reshape(HH5,100,100);

for i =1:100
pkw(i,ki) = kruskalwallis([H1(:,i),H2(:,i),H3(:,i),H4(:,i),H5(:,i)],{'bvic1','bvic2','aic','bic', 'aicc'},'off');
pan(i,ki) = anova1([H1(:,i),H2(:,i),H3(:,i),H4(:,i),H5(:,i)],{'bvic1','bvic2','aic','bic', 'aicc'},'off');
end
end

figure
hold on
boxplot(pkw);
ylabel('p-value','FontSize',16)
xlabel('forecast horizon','FontSize',16)
title('Kruskal-Wallis test','FontSize',16)
yline(0.05,'--k','LineWidth',1.3)
ylim([0 1.05])
xticks([0, 5, 10, 15, 20, 25, 30]) 
xticklabels([0, 5, 10, 15, 20, 25, 30])
figure
hold on
boxplot(pan)
title('One-way ANOVA test','FontSize',16)
ylabel('p-value','FontSize',16)
xlabel('forecast horizon','FontSize',16)
yline(0.05,'--k','LineWidth',1.3)
ylim([0 1.05])
xticks([0, 5, 10, 15, 20, 25, 30]) 
xticklabels([0, 5, 10, 15, 20, 25, 30])
