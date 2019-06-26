%% ---------------------- FILE INFO ---------------------------------------
%  Autor: Patryk Okoñski
%  Nazwa: mpc_porownanie.m
%  Przeznaczenie: Symulacja porównywanych algorytmów.
% -------------------------------------------------------------------------

clear all;
close all;
clc;

%% Parametry predykcyjne
% GPC
gpc_N = 10;
gpc_Nu = 3;
gpc_lambda_vec = [0.01];
% MPC-NO
no_N = 10;
no_Nu = 3;
no_lambda = 1;
% MPC-NPSL
npsl_N = no_N;
npsl_Nu = no_Nu;
npsl_lambda = no_lambda;
%% Alokacja wektorów
% GPC
gpc_y_zad_h = [];
gpc_y_zad_pH = [];
gpc_yp_h = [];
gpc_yp_pH = [];
gpc_q1 = [];
gpc_q3 = [];
gpc_e = [];
% MPC-NO
no_y_zad_h = [];
no_y_zad_pH = [];
no_yp_h = [];
no_yp_pH = [];
no_q1 = [];
no_q3 = [];
no_e = [];
% MPC-NPSL
npsl_y_zad_h = [];
npsl_y_zad_pH = [];
npsl_yp_h = [];
npsl_yp_pH = [];
npsl_q1 = [];
npsl_q3 = [];
npsl_e = [];
% Macierz b³êdów i czasów wykonania
met = zeros(10,max(size(gpc_lambda_vec)));

%% Wywo³anie algorytmów
% MPC-NO
disp('*Algorytm MPC-NO...');
start_no = tic;
[no_y_zad_h, no_y_zad_pH, no_yp_h, no_yp_pH, no_q1, no_q3, no_e] = c_mpc_no(no_N, no_Nu, no_lambda);
no_t = toc(start_no);
% MPC-NPSL
disp('*Algorytm MPC-NPSL...');
start_npsl = tic;
[npsl_y_zad_h, npsl_y_zad_pH, npsl_yp_h, npsl_yp_pH, npsl_q1, npsl_q3, npsl_e] = c_mpc_npsl(npsl_N, npsl_Nu, npsl_lambda);
npsl_t = toc(start_npsl);

%% Pêtla dla ró¿nych wartoœci gpc_lambda
for i = 1:max(size(gpc_lambda_vec))
    % GPC
    disp(['*Algorytm GPC, iteracja ',num2str(i),'...']);
    gpc_lambda = gpc_lambda_vec(i);
    start_gpc = tic;
    [gpc_y_zad_h, gpc_y_zad_pH, gpc_yp_h, gpc_yp_pH, gpc_q1, gpc_q3, gpc_e] = c_gpc(gpc_N, gpc_Nu, gpc_lambda);
    gpc_t = toc(start_gpc);
    
    % Inne
    cy = min([max(size(no_yp_h)) max(size(npsl_yp_h)) max(size(gpc_yp_h))]);
    cq = min([max(size(no_q1)) max(size(npsl_q1)) max(size(gpc_q1))]);
    y_zad_h = gpc_y_zad_h;
    y_zad_pH = gpc_y_zad_pH;
    
    % Zapisanie do 'met'
    met(1,i) = gpc_lambda_vec(i);
    met(2,i) = gpc_e(1);
    met(3,i) = gpc_e(2);
    met(4,i) = no_e(1);
    met(5,i) = no_e(2);
    met(6,i) = npsl_e(1);
    met(7,i) = npsl_e(2);
    met(8,i) = gpc_t;
    met(9,i) = no_t;
    met(10,i) = npsl_t;

    %% Wykresy
    % Wykres h
    w1 = figure(1);
    stairs(y_zad_h(1:cy),'r');
    hold on;
    stairs(no_yp_h(1:cy),'b');
    stairs(npsl_yp_h(1:cy),'g');
    stairs(gpc_yp_h(1:cy),'m');
    hold off;
    xlabel('Krok [k]');
    ylabel('h_{zad}, h [cm]');
    legend('h_{zad}','h_{NO}','h_{NPSL}','h_{GPC}');
    title(['MSE_{NO}= ',num2str(no_e(1)),',     MSE_{NPSL}= ',num2str(npsl_e(1)),',     MSE_{GPC}= ',num2str(gpc_e(1))]);
    grid on

    % Wykres pH
    w2 = figure(2);
    stairs(y_zad_pH(1:cy),'r');
    hold on;
    stairs(no_yp_pH(1:cy),'b');
    stairs(npsl_yp_pH(1:cy),'g');
    stairs(gpc_yp_pH(1:cy),'m');
    hold off;
    xlabel('Krok [k]');
    ylabel('pH_{zad}, pH [-]');
    legend('pH_{zad}','pH_{NO}','pH_{NPSL}','pH_{GPC}','Location','southeast');
    title(['MSE_{NO}= ',num2str(no_e(2)),',     MSE_{NPSL}= ',num2str(npsl_e(2)),',     MSE_{GPC}= ',num2str(gpc_e(2))]);
    grid on

    % Wykres q1
    w3 = figure(3);
    plot(no_q1(1:cq),'b');
    hold on
    plot(npsl_q1(1:cq),'g');
    plot(gpc_q1(1:cq),'m');
    hold off;
    xlabel('Krok [k]');
    ylabel('q1 [ml/s]');
    legend('q_1^{NO}','q_1^{NPSL}','q_1^{GPC}');
    title('q_1(k)');
    grid on

    % Wykres q3
    w4 = figure(4);
    plot(no_q3(1:cq),'b');
    hold on
    plot(npsl_q3(1:cq),'g');
    plot(gpc_q3(1:cq),'m');
    hold off;
    xlabel('Krok [k]');
    ylabel('q3 [ml/s]');
    legend('q_3^{NO}','q_3^{NPSL}','q_3^{GPC}','Location','southeast');
    title('q_3(k)');
    grid on


    %% Eksport wykresów
    % Budowanie nazwy pliku
%     w1_txt = sprintf('Wykresy/porownanie_h_N_%d_Nu_%d_lc_%0.2f_lgpc_%0.3f',no_N,no_Nu,no_lambda,gpc_lambda);
%     w2_txt = sprintf('Wykresy/porownanie_pH_N_%d_Nu_%d_lc_%0.2f_lgpc_%0.3f',no_N,no_Nu,no_lambda,gpc_lambda);
%     w3_txt = sprintf('Wykresy/porownanie_q1_N_%d_Nu_%d_lc_%0.2f_lgpc_%0.3f',no_N,no_Nu,no_lambda,gpc_lambda);
%     w4_txt = sprintf('Wykresy/porownanie_q3_N_%d_Nu_%d_lc_%0.2f_lgpc_%0.3f',no_N,no_Nu,no_lambda,gpc_lambda);
    
    % Budowanie nazwy pliku (Latex nie lubi kropek ...)
%     w1_txt = sprintf('Wykresy/porownanie_h_%d',i);
%     w2_txt = sprintf('Wykresy/porownanie_pH_%d',i);
%     w3_txt = sprintf('Wykresy/porownanie_q1_%d',i);
%     w4_txt = sprintf('Wykresy/porownanie_q3_%d',i);
    
    % Eksport
%     export_fig(w1, w1_txt, '-pdf', '-transparent', '-nocrop');
%     export_fig(w2, w2_txt, '-pdf', '-transparent', '-nocrop');
%     export_fig(w3, w3_txt, '-pdf', '-transparent', '-nocrop');
%     export_fig(w4, w4_txt, '-pdf', '-transparent', '-nocrop');
    
    % Czyszczenie wykresów dla kolejnej iteracji
%     clf(w1);
%     clf(w2);
%     clf(w3);
%     clf(w4);

end
% Eksport macierzy met
% save 'Dane/mpc_porownanie_met.mat' met;
% Zamkniêcie okien
% close(w1);
% close(w2);
% close(w3);
% close(w4);



