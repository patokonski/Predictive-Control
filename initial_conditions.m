%% ---------------------- FILE INFO ---------------------------------------
%  Autor: Patryk Okoñski
%  Nazwa: initial_conditions.m
%  Przeznaczenie: Generowanie warunków pocz¹tkowych.
% -------------------------------------------------------------------------

function [x] = initial_conditions(q1, q3)
    global Wa1 Wa2 Wa3 Wb1 Wb2 Wb3 K1 K2 A Cv q2 Ts
    
    Wa = (q1*Wa1 + q2*Wa2 + q3*Wa3)/(q1 + q2 + q3);
    Wb = (q1*Wb1 + q2*Wb2 + q3*Wb3)/(q1 + q2 + q3);
    h = ((q1 + q2 + q3)/Cv)^2;
    pH = calc_pH(Wa, Wb);
    x = [Wa Wb h pH];
end