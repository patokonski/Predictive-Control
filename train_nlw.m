%% ------------------------------------------------------------------------
%                       Praca dyplomowa magisterska
%  Autor: Patryk Okoñski
%  Temat: Nieliniowa regulacja predykcyjna procesu dwuwymiarowego
%  wykorzystuj¹ca neuronowy model Wienera.
%  ------------------------------------------------------------------------

%% ---------------------- FILE INFO ---------------------------------------
%  Nazwa: train_nlw.m
%  Przeznaczenie: Uczenie neuronoweogo modelu wienera dla wyjœæ procesu
%  neutralizacji.
% -------------------------------------------------------------------------

clear all;
clc;

init_reactor_params;

global Wa1 Wa2 Wa3 Wb1 Wb2 Wb3 K1 K2 A Cv q2 Ts

% Parametry
nb = 1;
nf = 1;
nk = 1;
K = 4;

% Parametry uczenia
nlhw_options = nlhwOptions;
nlhw_options.Display = 'on';
nlhw_options.SearchMethod = 'lm';
nlhw_options.SearchOption.MaxIter = 50;
nonlinearity_input = unitgain;
nonlinearity_output = sigmoidnet('NumberofUnits',K);
nonlinearity_output.LinearTerm= 'on';

% Rzêdy dynamiki
nb_vec = ones(1,2)*nb';
nf_vec = ones(1,2)*nf'; 
nk_vec = ones(1,2)*nk';  
n_vec = [nb_vec nf_vec nk_vec];

% Wczytanie danych
disp('Loading data ...');
load('Dane/train_data.mat');
load('Dane/validate_data.mat');
load('Dane/test_data.mat');

% Warunki pocz¹tkowe
q10 = 16.6;                             
q30 = 15.6;                             

x0 = initial_conditions(q10, q30);      
h0 = x0(3);
pH0 = x0(4);

%% Uczenie neuronowego modelu wienera
% Skalowanie danych
u1 = (train(:,1)-q10)/15;
u2 = (train(:,2)-q30)/15;
y_h = (train(:,3)-h0)/20;
y_pH = (train(:,4)-pH0)/5;

% Formatowanie
u = [u1 u2];
data_h = iddata(y_h, u, Ts);
data_pH = iddata(y_pH, u, Ts);

% Uczenie
model_h = nlhw(data_h, n_vec, nonlinearity_input, nonlinearity_output, nlhw_options);
model_pH = nlhw(data_pH, n_vec, nonlinearity_input, nonlinearity_output, nlhw_options);

% Symulacja
ymod_h = sim(model_h,  u);
ymod_pH = sim(model_pH, u);

% Skalowanie
ymod_h = ymod_h*20 + h0;
ymod_pH = ymod_pH*5 + pH0;

% Zapis
%save('Dane/model_h.mat','model_h');
%save('Dane/model_pH.mat','model_pH');


%% Weryfikacja
% Skalowanie danych
u1v = (validate(:,1)-q10)/15;
u2v = (validate(:,2)-q30)/15;

% Formatowanie
uv = [u1v u2v];

% Symulowanie
ymod_hv = sim(model_h, uv);
ymod_pHv = sim(model_pH, uv);

% Skalowanie
ymod_hv = ymod_hv*20 + h0;
ymod_pHv = ymod_pHv*5 + pH0;

%% Zapisanie
% save('Dane/nlhw_data.mat');


e1 = (train(:,3)-ymod_h)'*(train(:,3)-ymod_h);
e2 = (train(:,4)-ymod_pH)'*(train(:,4)-ymod_pH);
e3 = (validate(:,3)-ymod_hv)'*(validate(:,3)-ymod_hv);
e4 = (validate(:,4)-ymod_pHv)'*(validate(:,4)-ymod_pHv);

% e1
% e3
% e2
% e4

%% Wykresy K=3
% figure(1)
% plot(train(:,3));             
% hold on
% plot(ymod_h,'--r');                 
% hold off
% xlabel('Krok [k]');
% ylabel('h(k)');
% legend('proces', 'model');
% title(['Uczenie wyjœcia h, SSE = ',num2str(e1)]);
% export_fig uczenie_nlhw_h_K3.pdf -transparent -nocrop
% 
% figure(2)
% plot(train(:,4));             
% hold on
% plot(ymod_pH,'--r');                
% hold off
% xlabel('Krok [k]');
% ylabel('pH(k)');
% legend('proces', 'model');
% title(['Uczenie wyjœcia pH, SSE = ',num2str(e2)]);
% export_fig uczenie_nlhw_pH_K3.pdf -transparent -nocrop
% 
% figure(3)
% plot(validate(:,3));             
% hold on
% plot(ymod_hv,'--r');                 
% hold off
% xlabel('Krok [k]');
% ylabel('h(k)');
% legend('proces', 'model');
% title(['Weryfikacja wyjœcia h, SSE = ',num2str(e3)]);
% export_fig weryfikacja_nlhw_h_K3.pdf -transparent -nocrop
% 
% figure(4)
% plot(validate(:,4));             
% hold on
% plot(ymod_pHv,'--r');                
% hold off
% xlabel('Krok [k]');
% ylabel('pH(k)');
% legend('proces', 'model');
% title(['Weryfikacja wyjœcia pH, SSE = ',num2str(e4)]);
% export_fig weryfikacja_nlhw_pH_K3.pdf -transparent -nocrop

%% Wykresy K=5
% figure(1)
% plot(train(:,3));             
% hold on
% plot(ymod_h,'--r');                 
% hold off
% xlabel('Krok [k]');
% ylabel('h(k)');
% legend('proces', 'model');
% title(['Uczenie wyjœcia h, SSE = ',num2str(e1)]);
% export_fig uczenie_nlhw_h_K5.pdf -transparent -nocrop
% 
% figure(2)
% plot(train(:,4));             
% hold on
% plot(ymod_pH,'--r');                
% hold off
% xlabel('Krok [k]');
% ylabel('pH(k)');
% legend('proces', 'model');
% title(['Uczenie wyjœcia pH, SSE = ',num2str(e2)]);
% export_fig uczenie_nlhw_pH_K5.pdf -transparent -nocrop
% 
% figure(3)
% plot(validate(:,3));             
% hold on
% plot(ymod_hv,'--r');                 
% hold off
% xlabel('Krok [k]');
% ylabel('h(k)');
% legend('proces', 'model');
% title(['Weryfikacja wyjœcia h, SSE = ',num2str(e3)]);
% export_fig weryfikacja_nlhw_h_K5.pdf -transparent -nocrop
% 
% figure(4)
% plot(validate(:,4));             
% hold on
% plot(ymod_pHv,'--r');                
% hold off
% xlabel('Krok [k]');
% ylabel('pH(k)');
% legend('proces', 'model');
% title(['Weryfikacja wyjœcia pH, SSE = ',num2str(e4)]);
% export_fig weryfikacja_nlhw_pH_K5.pdf -transparent -nocrop

%% Wykresy K=10
% figure(1)
% plot(train(:,3));             
% hold on
% plot(ymod_h,'--r');                 
% hold off
% xlabel('Krok [k]');
% ylabel('h(k)');
% legend('proces', 'model');
% title(['Uczenie wyjœcia h, SSE = ',num2str(e1)]);
% export_fig uczenie_nlhw_h_K10.pdf -transparent -nocrop
% 
% figure(2)
% plot(train(:,4));             
% hold on
% plot(ymod_pH,'--r');                
% hold off
% xlabel('Krok [k]');
% ylabel('pH(k)');
% legend('proces', 'model');
% title(['Uczenie wyjœcia pH, SSE = ',num2str(e2)]);
% export_fig uczenie_nlhw_pH_K10.pdf -transparent -nocrop
% 
% figure(3)
% plot(validate(:,3));             
% hold on
% plot(ymod_hv,'--r');                 
% hold off
% xlabel('Krok [k]');
% ylabel('h(k)');
% legend('proces', 'model');
% title(['Weryfikacja wyjœcia h, SSE = ',num2str(e3)]);
% export_fig weryfikacja_nlhw_h_K10.pdf -transparent -nocrop
% 
% figure(4)
% plot(validate(:,4));             
% hold on
% plot(ymod_pHv,'--r');                
% hold off
% xlabel('Krok [k]');
% ylabel('pH(k)');
% legend('proces', 'model');
% title(['Weryfikacja wyjœcia pH, SSE = ',num2str(e4)]);
% export_fig weryfikacja_nlhw_pH_K10.pdf -transparent -nocrop

