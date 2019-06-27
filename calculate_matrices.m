%% ---------------------- FILE INFO ---------------------------------------
%  Autor: Patryk Okoński
%  Nazwa: calculate_matrices.m
%  Przeznaczenie: Obliczenie macierzy S, G i K na potrzeby algorytmu GPC.
% -------------------------------------------------------------------------
function [G, K] = calculate_matrices(Hp,Hs,wsplambda,wspmi,ny,nu,na,nb)
    %% MODEL LINIOWY
     load model_linear.mat;
%     load model_linear_nlhw.mat;
    %% Alokacja macierzy
    s = zeros(2,2,Hp);

    S = [];

    G = [];

    K = [];

    %% Obliczenie wszystkich macierzy S wsp odpowiedzi skokowej s
    % Format S(:,:,j) = [s(1,1,j) s(1,2,j); s(2,1,j) s(2,2,j)]   s(m,n,j) = wartosc;
    for j = 1:Hp
        for lm = 1:ny
            for ln = 1:nu
              for i = 1:min(j,nb)
                  s(lm,ln,j) = s(lm,ln,j) + b(lm,ln,i);
              end
              for i = 1:min(j-1,na)
                  s(lm,ln,j) = s(lm,ln,j) - a(lm,i)*s(lm,ln,j-i);
              end
            end
        end
        S(:,:,j) = [s(1,1,j) s(1,2,j); s(2,1,j) s(2,2,j)];
    end
    Sgpc = S;
	save SGPC.mat Sgpc 	
    % Obliczenie macierzy G
    for i = 1:Hs
        for j = i:Hp
            G((2*j-1):2*j,(2*i-1):2*i) = S(:,:,j-i+1);
        end
    end
   
    %% Macierze M i L
    i=1;
    for p=1:Hp
        for m=1:2
            MI(i,i)=wspmi(m); i=i+1;
        end
    end
    i=1;
    for p=1:Hs
        for n=1:2
            LAMBDA(i,i)=wsplambda(n); i=i+1;
        end
    end

    %% wzmocnienie K
    K=inv(G'*MI*G+LAMBDA)*G'*MI;
    % print S
end
