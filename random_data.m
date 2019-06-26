%% ------------------------------------------------------------------------
%                       Praca dyplomowa magisterska
%  Autor: Patryk Oko�ski
%  Temat: Nieliniowa regulacja predykcyjna procesu dwuwymiarowego
%  wykorzystuj�ca neuronowy model Wienera.
%  ------------------------------------------------------------------------

%% ---------------------- FILE INFO ---------------------------------------
%  Nazwa: random_data.m
%  Przeznaczenie: Generowanie losowych sterowa�.
% -------------------------------------------------------------------------

function [data] = random_data(q1, q3, q1min, q1max, q3min, q3max, qmod, iter)
    data(1:qmod,1) = q1;
    data(1:qmod,2) = q3;
    for i = (qmod+1):iter
        if (mod(i,qmod) == 0)
            q1 = q1min + (q1max - q1min)*rand(1);
            q3 = q3min + (q3max - q3min)*rand(1);
        end
        data(i,1) = q1;
        data(i,2) = q3;
    end
end