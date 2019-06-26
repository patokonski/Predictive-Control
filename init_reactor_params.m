%% ---------------------- FILE INFO ---------------------------------------
%  Autor: Patryk Okoñski
%  Nazwa: init_reactor_params.m
%  Przeznaczenie: Inicjalizacja parametrów reaktora.
% -------------------------------------------------------------------------

function init_reactor_params
    global Wa1 Wa2 Wa3 Wb1 Wb2 Wb3 K1 K2 A Cv q2 Ts
    
    Wa1 = 0.003;                %[M]
    Wa2 = -0.03;                %[M]
    Wa3 = -0.00305;             %[M]
    Wb1 = 0;                    %[M]          
    Wb2 = 0.03;                 %[M]
    Wb3 = 0.00005;              %[M]
    K1 = -log10(4.47e-7);       %[-]
    K2 = -log10(5.62e-11);      %[-]
    A = 207;                    %[cm2]
    Cv = 8.75;                  %[ml/cm s]
    q2 = 0.55;                  %[ml/s]
    Ts = 10;                    %[s]
end