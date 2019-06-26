%% ------------------------------------------------------------------------
%                       Praca dyplomowa magisterska
%  Autor: Patryk Okoñski
%  Temat: Nieliniowa regulacja predykcyjna procesu dwuwymiarowego
%  wykorzystuj¹ca neuronowy model Wienera.
%  ------------------------------------------------------------------------

%% ---------------------- FILE INFO ---------------------------------------
%  Nazwa: calc_pH.m
%  Przeznaczenie: Liczenie wartoœci pH.
% -------------------------------------------------------------------------

function [pH] = calc_pH(Wa, Wb)
    global Wa1 Wa2 Wa3 Wb1 Wb2 Wb3 K1 K2 A Cv q2 Ts
    
    p4 = 10^(-14)*10^(-K2);
    p3 = ( Wa*10^(-K2)+10^(-14)+Wb*2*10^(-K2) );
    p2 = ( Wa+10^(-14)*10^K1+Wb-10^(-K2) );
    p1 = ( Wa*10^(K1)-1 );
    p0 = -10^K1;
    pH = log10(max(roots([p4 p3 p2 p1 p0])));
    
end