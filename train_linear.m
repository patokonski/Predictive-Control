%% ---------------------- FILE INFO ---------------------------------------
%  Autor: Patryk Okoñski
%  Nazwa: train_linear.m
%  Przeznaczenie: Wyznaczenie modelu liniowego procesu neutralizacji.
% -------------------------------------------------------------------------

clear all;
clc;

global Wa1 Wa2 Wa3 Wb1 Wb2 Wb3 K1 K2 A Cv q2 Ts

init_reactor_params;
load('Dane\train_data.mat');

q10 = 16.6;
q30 = 15.6;

x0 = initial_conditions(q10, q30);
Wa0 = x0(1);
Wb0 = x0(2);
h0 = x0(3);
pH0 = x0(4);

u1 = (train(:,1)-q10)/15;
u2 = (train(:,2)-q30)/15;
yh = (train(:,3)-h0)/20;
ypH = (train(:,4)-pH0)/5;

kstart = 3;
kend = length(u1);
ny = 2;
nu = 2;
na = 2;
nb = 2;

Mh = [u1(kstart-1:kend-1) u1(kstart-2:kend-2) u2(kstart-1:kend-1) u2(kstart-2:kend-2) yh(kstart-1:kend-1) yh(kstart-2:kend-2)];
MpH = [u1(kstart-1:kend-1) u1(kstart-2:kend-2) u2(kstart-1:kend-1) u2(kstart-2:kend-2) ypH(kstart-1:kend-1) ypH(kstart-2:kend-2)];   

wh = Mh\yh(kstart:kend);
wpH = MpH\ypH(kstart:kend);

b(1,1,1) = wh(1);
b(1,1,2) = wh(2);
b(1,2,1) = wh(3);
b(1,2,2) = wh(4);
b(2,1,1) = wpH(1);
b(2,1,2) = wpH(2);
b(2,2,1) = wpH(3);
b(2,2,2) = wpH(4);
a(1,1) = -wh(5);
a(1,2) = -wh(6);
a(2,1) = -wpH(5);
a(2,2) = -wpH(6);    

u(1,1:length(u1)) = u1;
u(2,1:length(u1)) = u2;

y(1,1:kstart-1) = yh(1:kstart-1);
y(2,1:kstart-1) = ypH(1:kstart-1);

for k = kstart:kend
    for m=1:ny
        y(m,k)=0.0;
        for n=1:nu
            for i=1:nb
                y(m,k)=y(m,k)+b(m,n,i)*u(n,k-i);
            end
        end
        for i=1:na
            y(m,k)=y(m,k)-a(m,i)*y(m,k-i);
        end
    end
end

ymodh = 20*y(1,1:end) + h0;
ymodpH = 5*y(2,1:end) + pH0;

save Dane\model_linear a b

figure(1)
plot(train(:,3),'b')
xlabel('k')
ylabel('h')
hold on
plot(ymodh, '--r')
title('h, ymod')
legend('h','ymod');

figure(2)
plot(train(:,4),'b')
xlabel('k')
ylabel('pH')
hold on
plot(ymodpH, '--r')    
title('pH, ymod')
legend('pH','ymod')
