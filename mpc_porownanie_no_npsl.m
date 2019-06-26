%% ---------------------- FILE INFO ---------------------------------------
%  Autor: Patryk Okoñski
%  Nazwa: mpc_porownanie.m
%  Przeznaczenie: Symulacja porównywanych algorytmów.
% -------------------------------------------------------------------------

clear all;
close all;
clc;

run kolory;

%% Parametry predykcyjne
% MPC-NO
no_N = 10;
no_Nu = 3;
% no_lambda = 1;
% MPC-NPSL
npsl_N = no_N;
npsl_Nu = no_Nu;
% npsl_lambda = no_lambda;
% wektor lambda
vlam = [0.1 0.25 0.5 1];
%% Alokacja wektorów

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
% Macierz b³êdów i czasów wykonania (14,x,max(size(vlam)))
met = [];
mt = [];
%% Pêtla dla ró¿nych wartoœci lambda
for i = 1:max(size(vlam))
    no_lambda = vlam(i)
    npsl_lambda = no_lambda;
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
    
    mt(1,i) = no_t;
    mt(2,i) = npsl_t;
    
    % Inne
    cy = min([max(size(no_yp_h)) max(size(npsl_yp_h))]);
    cq = min([max(size(no_q1)) max(size(npsl_q1))]);
    y_zad_h = no_y_zad_h;
    y_zad_pH = no_y_zad_pH;
    
    % Zapisanie do 'met'
%     met(1,:,i) = y_zad_h(1:cy);
%     met(2,:,i) = y_zad_pH(1:cy);
%     met(3,:,i) = no_yp_h(1:cy);
%     met(4,:,i) = no_yp_pH(1:cy);
%     met(5,:,i) = no_q1;
%     met(6,:,i) = no_q3;
%     met(7,:,i) = no_e(1);
%     met(8,:,i) = no_e(2);
%     met(9,:,i) = npsl_yp_h;
%     met(10,:,i) = npsl_yp_pH;
%     met(11,:,i) = npsl_q1;
%     met(12,:,i) = npsl_q3;
%     met(13,:,i) = npsl_e(1);
%     met(14,:,i) = npsl_e(2)
%     %% Wykresy
%     % Wykres h
%     w1 = figure(1);
%     stairs(y_zad_h(1:cy),'r');
%     hold on;
%     stairs(no_yp_h(1:cy),styl1,'color',kolor1,'linewidth',grubosc);
%     stairs(npsl_yp_h(1:cy),styl3,'color',kolor3,'linewidth',grubosc);
%     hold off;
%     xlabel('Krok [k]');
%     ylabel('h_{zad}, h [cm]');
%     legend('h_{zad}','h_{NO}','h_{NPSL}');
%     title(['SSE_{NO}= ',num2str(no_e(1)),',     SSE_{NPSL}= ',num2str(npsl_e(1)),', N = ',num2str(no_N),', N_u = ',num2str(no_Nu),', \lambda = ',num2str(vlam(i))]);
%     grid on
% 
%     % Wykres pH
%     w2 = figure(2);
%     stairs(y_zad_pH(1:cy),'r');
%     hold on;
%     stairs(no_yp_pH(1:cy),styl1,'color',kolor1,'linewidth',grubosc);
%     stairs(npsl_yp_pH(1:cy),styl3,'color',kolor3,'linewidth',grubosc);
%     hold off;
%     xlabel('Krok [k]');
%     ylabel('pH_{zad}, pH [-]');
%     legend('pH_{zad}','pH_{NO}','pH_{NPSL}','Location','southeast');
%     title(['SSE_{NO}= ',num2str(no_e(2)),',     SSE_{NPSL}= ',num2str(npsl_e(2)),', N = ',num2str(no_N),', N_u = ',num2str(no_Nu),', \lambda = ',num2str(vlam(i))]);
%     grid on
% 
%     % Wykres q1
%     w3 = figure(3);
%     plot(no_q1(1:cq),styl2,'color',kolor2,'linewidth',grubosc);
%     hold on
%     plot(npsl_q1(1:cq),styl3,'color',kolor3,'linewidth',grubosc);
%     hold off;
%     xlabel('Krok [k]');
%     ylabel('q1 [ml/s]');
%     legend('q_1^{NO}','q_1^{NPSL}');
%     title(['q_1(k), N = ',num2str(no_N),', N_u = ',num2str(no_Nu),', \lambda = ',num2str(vlam(i))]);
%     grid on
% 
%     % Wykres q3
%     w4 = figure(4);
%     plot(no_q3(1:cq),styl2,'color',kolor2,'linewidth',grubosc);
%     hold on
%     plot(npsl_q3(1:cq),styl3,'color',kolor3,'linewidth',grubosc);
%     hold off;
%     xlabel('Krok [k]');
%     ylabel('q3 [ml/s]');
%     legend('q_3^{NO}','q_3^{NPSL}','Location','southeast');
%     title(['q_3(k), N = ',num2str(no_N),', N_u = ',num2str(no_Nu),', \lambda = ',num2str(vlam(i))]);
%     grid on
%     
%     clearvars no_y_zad_h no_y_zad_pH no_yp_h no_yp_pH no_q1 no_q3 no_e
%     clearvars npsl_y_zad_h npsl_y_zad_pH npsl_yp_h npsl_yp_pH npsl_q1 npsl_q3 npsl_e
    
%     %% Eksport wykresów
%     % Budowanie nazwy pliku
% %     w1_txt = sprintf('Wykresy/porownanie_h_N_%d_Nu_%d_lc_%0.2f',no_N,no_Nu,no_lambda);
% %     w2_txt = sprintf('Wykresy/porownanie_pH_N_%d_Nu_%d_lc_%0.2f',no_N,no_Nu,no_lambda);
% %     w3_txt = sprintf('Wykresy/porownanie_q1_N_%d_Nu_%d_lc_%0.2f',no_N,no_Nu,no_lambda);
% %     w4_txt = sprintf('Wykresy/porownanie_q3_N_%d_Nu_%d_lc_%0.2f',no_N,no_Nu,no_lambda);
%     
%     % Budowanie nazwy pliku (Latex nie lubi kropek ...)
%     w1_txt = sprintf('Wykresy/porownanie_h_%d',i);
%     w2_txt = sprintf('Wykresy/porownanie_pH_%d',i);
%     w3_txt = sprintf('Wykresy/porownanie_q1_%d',i);
%     w4_txt = sprintf('Wykresy/porownanie_q3_%d',i);
%     
%     % Eksport
%     export_fig(w1, w1_txt, '-pdf', '-transparent', '-nocrop');
%     export_fig(w2, w2_txt, '-pdf', '-transparent', '-nocrop');
%     export_fig(w3, w3_txt, '-pdf', '-transparent', '-nocrop');
%     export_fig(w4, w4_txt, '-pdf', '-transparent', '-nocrop');
%     export_fig(w1, w1_txt, '-png', '-transparent', '-nocrop');
%     export_fig(w2, w2_txt, '-png', '-transparent', '-nocrop');
%     export_fig(w3, w3_txt, '-png', '-transparent', '-nocrop');
%     export_fig(w4, w4_txt, '-png', '-transparent', '-nocrop');
%     
%     % Czyszczenie wykresów dla kolejnej iteracji
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



