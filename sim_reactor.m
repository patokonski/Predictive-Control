%% ---------------------- FILE INFO ---------------------------------------
%  Autor: Patryk Okoński
%  Nazwa: sim_reactor.m
%  Przeznaczenie: Rozwiązanie równań różniczkowych procesu z ode45.
% -------------------------------------------------------------------------

function xd = sim_reactor(t, x, q1, q3)
    global Wa1 Wa2 Wa3 Wb1 Wb2 Wb3 K1 K2 A Cv q2 Ts
    xd = [0; 0; 0];
    Wa = x(1);
    Wb = x(2);
    h = x(3);
    xd(1) = (q1*Wa1 + q2*Wa2 + q3*Wa3 - Wa*(q1 + q2 + q3))/(A*h);
    xd(2) = (q1*Wb1 + q2*Wb2 + q3*Wb3 - Wb*(q1 + q2 + q3))/(A*h);
    xd(3) = (q1 + q2 + q3 - Cv*sqrt(h))/A;
end
