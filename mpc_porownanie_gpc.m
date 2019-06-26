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
gpc_lambda_vec = [0.001 0.01 0.1 1];
%% Alokacja wektorów
% GPC
gpc_y_zad_h = [];
gpc_y_zad_pH = [];
gpc_yp_h = [];
gpc_yp_pH = [];
gpc_q1 = [];
gpc_q3 = [];
gpc_e = [];
% Macierz b³êdów i czasów wykonania
mgpcy = [];
mgpce = [];

%% Pêtla dla ró¿nych wartoœci gpc_lambda
for i = 1:max(size(gpc_lambda_vec))   
    % GPC
    disp(['*Algorytm GPC, iteracja ',num2str(i),'...']);
    gpc_lambda = gpc_lambda_vec(i);
    start_gpc = tic;
    [gpc_y_zad_h, gpc_y_zad_pH, gpc_yp_h, gpc_yp_pH, gpc_q1, gpc_q3, gpc_e] = c_gpc(gpc_N, gpc_Nu, gpc_lambda);
    gpc_t = toc(start_gpc);
    
    % Inne
%     cy = min([max(size(no_yp_h)) max(size(npsl_yp_h)) max(size(gpc_yp_h))]);
%     cq = min([max(size(no_q1)) max(size(npsl_q1)) max(size(gpc_q1))]);
    if (i == 1)
        y_zad_h = gpc_y_zad_h;
        y_zad_pH = gpc_y_zad_pH;
    end
    
    mgpcy(1,:,i) = gpc_yp_h;
    mgpcy(2,:,i) = gpc_yp_pH;
    mgpcy(3,:,i) = gpc_q1;
    mgpcy(4,:,i) = gpc_q3;
    
    mgpce(1,i) = gpc_e(1);
    mgpce(2,i) = gpc_e(2);
    % Zapisanie do 'met'
%     met(1,i) = gpc_lambda_vec(i);
%     met(2,i) = gpc_e(1);
%     met(3,i) = gpc_e(2);
%     met(4,i) = no_e(1);
%     met(5,i) = no_e(2);
%     met(6,i) = npsl_e(1);
%     met(7,i) = npsl_e(2);
%     met(8,i) = gpc_t;
%     met(9,i) = no_t;
%     met(10,i) = npsl_t;
end

%% Wykresy
run kolory
% Wykres h
w1 = figure(1);
stairs(y_zad_h,'r');
hold on;
stairs(mgpcy(1,:,1),styl1,'color',kolor1,'linewidth',grubosc);
stairs(mgpcy(1,:,2),styl2,'color',kolor2,'linewidth',grubosc);
stairs(mgpcy(1,:,3),styl3,'color',kolor3,'linewidth',grubosc);
stairs(mgpcy(1,:,4),styl4,'color',kolor4,'linewidth',grubosc);
hold off;
xlabel('Krok [k]');
ylabel('h_{zad}, h [cm]');
legend('h_{zad}','\lambda=0.001','\lambda=0.01','\lambda=0.1','\lambda=1');
title(['SSE_{\lambda=0.001}= ',num2str(mgpce(1,1)),'; SSE_{\lambda=0.01}= ',num2str(mgpce(1,2)),'; SSE_{\lambda=0.1}= ',num2str(mgpce(1,3)),'; SSE_{\lambda=1}= ',num2str(mgpce(1,4))]);
grid on

% Wykres pH
w2 = figure(2);
stairs(y_zad_pH,'r');
hold on;
stairs(mgpcy(2,:,1),styl1,'color',kolor1,'linewidth',grubosc);
stairs(mgpcy(2,:,2),styl2,'color',kolor2,'linewidth',grubosc);
stairs(mgpcy(2,:,3),styl3,'color',kolor3,'linewidth',grubosc);
stairs(mgpcy(2,:,4),styl4,'color',kolor4,'linewidth',grubosc);
hold off;
xlabel('Krok [k]');
ylabel('pH_{zad}, pH [-]');
legend('pH_{zad}','\lambda=0.001','\lambda=0.01','\lambda=0.1','\lambda=1','Location','southeast');
title(['SSE_{\lambda=0.001}= ',num2str(mgpce(2,1)),'; SSE_{\lambda=0.01}= ',num2str(mgpce(2,2)),'; SSE_{\lambda=0.1}= ',num2str(mgpce(2,3)),'; SSE_{\lambda=1}= ',num2str(mgpce(2,4))]);
grid on

% Wykres q1
w3 = figure(3);
stairs(mgpcy(3,:,1),styl1,'color',kolor1,'linewidth',grubosc);
hold on
stairs(mgpcy(3,:,2),styl2,'color',kolor2,'linewidth',grubosc);
stairs(mgpcy(3,:,3),styl3,'color',kolor3,'linewidth',grubosc);
stairs(mgpcy(3,:,4),styl4,'color',kolor4,'linewidth',grubosc);
hold off;
xlabel('Krok [k]');
ylabel('q1 [ml/s]');
legend('\lambda=0.001',' \lambda=0.01',' \lambda=0.1',' \lambda=1');
title('q_1(k)');
grid on

% Wykres q3
w4 = figure(4);
stairs(mgpcy(4,:,1),styl1,'color',kolor1,'linewidth',grubosc);
hold on
stairs(mgpcy(4,:,2),styl2,'color',kolor2,'linewidth',grubosc);
stairs(mgpcy(4,:,3),styl3,'color',kolor3,'linewidth',grubosc);
stairs(mgpcy(4,:,4),styl4,'color',kolor4,'linewidth',grubosc);
hold off;
xlabel('Krok [k]');
ylabel('q3 [ml/s]');
legend('\lambda=0.001',' \lambda=0.01',' \lambda=0.1',' \lambda=1');
title('q_3(k)');
grid on


%% Eksport wykresów
% Budowanie nazwy pliku
% w1_txt = sprintf('Wykresy/porownanie_h_N_%d_Nu_%d_lc_%0.2f_lgpc_%0.3f',no_N,no_Nu,no_lambda,gpc_lambda);
% w2_txt = sprintf('Wykresy/porownanie_pH_N_%d_Nu_%d_lc_%0.2f_lgpc_%0.3f',no_N,no_Nu,no_lambda,gpc_lambda);
% w3_txt = sprintf('Wykresy/porownanie_q1_N_%d_Nu_%d_lc_%0.2f_lgpc_%0.3f',no_N,no_Nu,no_lambda,gpc_lambda);
% w4_txt = sprintf('Wykresy/porownanie_q3_N_%d_Nu_%d_lc_%0.2f_lgpc_%0.3f',no_N,no_Nu,no_lambda,gpc_lambda);

% Budowanie nazwy pliku (Latex nie lubi kropek ...)
w1_txt = sprintf('Wykresy/gpc_h');
w2_txt = sprintf('Wykresy/gpc_pH');
w3_txt = sprintf('Wykresy/gpc_q1');
w4_txt = sprintf('Wykresy/gpc_q3');

% Eksport
% export_fig(w1, w1_txt, '-pdf', '-transparent', '-nocrop');
% export_fig(w2, w2_txt, '-pdf', '-transparent', '-nocrop');
% export_fig(w3, w3_txt, '-pdf', '-transparent', '-nocrop');
% export_fig(w4, w4_txt, '-pdf', '-transparent', '-nocrop');
% export_fig(w1, w1_txt, '-png', '-transparent', '-nocrop');
% export_fig(w2, w2_txt, '-png', '-transparent', '-nocrop');
% export_fig(w3, w3_txt, '-png', '-transparent', '-nocrop');
% export_fig(w4, w4_txt, '-png', '-transparent', '-nocrop');

% Czyszczenie wykresów dla kolejnej iteracji
% clf(w1);
% clf(w2);
% clf(w3);
% clf(w4);
    
    
% Eksport macierzy met
% save 'Dane/mpc_porownanie_met.mat' met;
% Zamkniêcie okien
% close(w1);
% close(w2);
% close(w3);
% close(w4);



