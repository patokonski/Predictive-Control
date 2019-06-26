%% ---------------------- FILE INFO ---------------------------------------
%  Autor: Patryk Okoñski
%  Nazwa: generate_h_pH.m
%  Przeznaczenie: Generowanie wyjœcia h i pH dla danych sterowañ.
% -------------------------------------------------------------------------

function [h_pH] = generate_h_pH(data, x0)
    global Wa1 Wa2 Wa3 Wb1 Wb2 Wb3 K1 K2 A Cv q2 Ts
    q1 = data(:,1);
    q3 = data(:,2);
    for k = 2:length(q1)
        [t,x] = ode45(@(t,x) sim_reactor(t, x, q1(k-1), q3(k-1)), [(k-1)*Ts (k)*Ts],x0);
        x0 = [x(end,1) x(end,2) x(end,3)];
        h(k) = x(end,3);
        pH(k) = calc_pH(x(end,1), x(end,2));
    end
    h_pH(:,1) = h;
    h_pH(:,2) = pH;
end