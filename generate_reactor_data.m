%% ------------------------------------------------------------------------
%                       Praca dyplomowa magisterska
%  Autor: Patryk Okoñski
%  Temat: Nieliniowa regulacja predykcyjna procesu dwuwymiarowego
%  wykorzystuj¹ca neuronowy model Wienera.
%  ------------------------------------------------------------------------

%% ---------------------- FILE INFO ---------------------------------------
%  Nazwa: generate_reactor_data.m
%  Przeznaczenie: Generowanie zbiorów danych z reaktora.
% -------------------------------------------------------------------------

global Wa1 Wa2 Wa3 Wb1 Wb2 Wb3 K1 K2 A Cv q2 Ts

clear all;
clc;

init_reactor_params;                    

q10 = 16.6;                             
q30 = 15.6;                             

x0 = initial_conditions(q10, q30);      

q1min = 0;                              
q1max = 30;
q3min = q1min;
q3max = q1max;
qmod = 30;                              

train = random_data(q10, q30, q1min, q1max, q3min, q3max, qmod, 3000);
validate = random_data(q10, q30, q1min, q1max, q3min, q3max, qmod, 3000);
test = random_data(q10, q30, q1min, q1max, q3min, q3max, qmod, 1000);

train(:,3:4) = generate_h_pH(train, x0(1:3));
validate(:,3:4) = generate_h_pH(validate, x0(1:3));
test(:,3:4) = generate_h_pH(test, x0(1:3));

train(1,3:4) = x0(3:4);
validate(1,3:4) = x0(3:4);
test(1,3:4) = x0(3:4);

% save('Dane/train_data.mat','train');
% save('Dane/validate_data.mat','validate');
% save('Dane/test_data.mat','test');

figure(1)
plot(train(:,1))
grid on
xlabel('k');
ylabel('q_1');
% export_fig q1_czyste.pdf -transparent -nocrop

figure(2)
plot(train(:,2))
grid on
xlabel('k');
ylabel('q_3');
% export_fig q3_czyste.pdf -transparent -nocrop



figure(3)
plot(train(:,3))
grid on
xlabel('k');
ylabel('h(k)');
% export_fig h_czyste.pdf -transparent -nocrop

figure(4)
plot(train(:,4))   
grid on
xlabel('k');
ylabel('pH(k)');
% export_fig pH_czyste.pdf -transparent -nocrop