clear all 
close all
clc
rng default

N = 100;
alpha = 0;
 pol1 = [0.9; 0.6 + 0.6i; 0.6 - 0.6i; -0.1 + 0.35i; -0.1- 0.35i];
 pol2 = [0.2+0.2i; 0.2-0.2i; -0.2; -0.6 + 0.6i; -0.6 - 0.6i];
 pol3 = [-0.4+0.2i; -0.4-0.2i; 0.5; 0.1-0.9i; 0.1+0.9i];
 pol4 = [-0.3; -0.1; 0.3; 0.5 + 0.1i; 0.5 - 0.1i];



p1 = poly(pol1);
p2 = poly(pol2);
p3 = poly(pol3);
p4 = poly(pol4);


Nts = 300*N;

W = normrnd(0,1,[1,Nts]);
% X = simulate(mdl0,Nts);
X1 = filter(1,p1,W);
% X1 = (X1-mean(X1))/std(X1);
X2 = filter(1,p2,W);
X3 = filter(1,p3,W);
X4 = filter(1,p4,W);


z = tf('z');

s1 = tf(1,p1,-1,'Variable','z^-1');
s2 = tf(1,p2,-1,'Variable','z^-1');
s3 = tf(1,p3,-1,'Variable','z^-1');
s4 = tf(1,p4,-1,'Variable','z^-1');



[~,idx] = max(abs(pol2));
w0 = angle(pol2(idx))
w0deg = w0*180/pi

%% plots
% f = figure;
% zplane([],pol,'r')
% f.Position = [10 300 500 400];
% g = figure;
% impz(1,par,40)
% g.Position = [510 300 500 400];
% q = figure;
% plot(X)
% q.Position = [1010 300 500 400];

msz = 0.1;
msp = 12;
% figure
% subplot(1,4,1), pzplot(s1,'r')
% a = findobj(gca,'type','line');
% set(a(2),'markersize',msz)
% set(a(3),'markersize',msp)
% title('Case 1')
% xlabel('')
% subplot(1,4,2), pzplot(s2,'r')
% a = findobj(gca,'type','line');
% set(a(2),'markersize',msz)
% set(a(3),'markersize',msp)
% title('Case 2')
% xlabel('')
% ylabel('')
% subplot(1,4,3), pzplot(s3,'r')
% a = findobj(gca,'type','line');
% set(a(2),'markersize',msz)
% set(a(3),'markersize',msp)
% title('Case 3')
% subplot(1,4,4), pzplot(s4,'r')
% a = findobj(gca,'type','line');
% set(a(2),'markersize',msz)
% set(a(3),'markersize',msp)
% title('Case 4')
% ylabel('')
wind = 1576:1625;
ax = [-5 5];

figure
subplot(2,4,1), pzplot(s1,'r')
a = findobj(gca,'type','line');
set(a(2),'markersize',msz)
set(a(3),'markersize',msp)
set(a(3),'LineWidth',2)
title('Case 1')

subplot(2,4,2), pzplot(s2,'r')
a = findobj(gca,'type','line');
set(a(2),'markersize',msz)
set(a(3),'markersize',msp)
set(a(3),'LineWidth',2)
title('Case 2')
ylabel('')

subplot(2,4,3), pzplot(s3,'r')
a = findobj(gca,'type','line');
set(a(2),'markersize',msz)
set(a(3),'markersize',msp)
set(a(3),'LineWidth',2)
title('Case 3')
ylabel('')
subplot(2,4,4), pzplot(s4,'r')
a = findobj(gca,'type','line');
set(a(2),'markersize',msz)
set(a(3),'markersize',msp)
set(a(3),'LineWidth',2)
title('Case 4')
ylabel('')

subplot(2,4,5), plot(wind,X1(wind),'LineWidth',1.2)
ylim(ax*2)
ylabel('Signal Y_t')
xlabel('t')
xticklabels([])
subplot(2,4,6), plot(wind,X2(wind),'LineWidth',1.2)
ylim(ax)
xlabel('t')
xticklabels([])
subplot(2,4,7), plot(wind,X3(wind),'LineWidth',1.2)
ylim(ax)
xlabel('t')
xticklabels([])
subplot(2,4,8), plot(wind,X4(wind),'LineWidth',1.2)
ylim(ax)
xlabel('t')
xticklabels([])