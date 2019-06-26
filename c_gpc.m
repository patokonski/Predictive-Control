%% ------------------------------------------------------------------------
%                       Praca dyplomowa magisterska
%  Autor: Patryk Okoñski
%  Temat: Nieliniowa regulacja predykcyjna procesu dwuwymiarowego
%  wykorzystuj¹ca neuronowy model Wienera
%  ------------------------------------------------------------------------

%% ---------------------- FILE INFO ---------------------------------------
%  Nazwa: c_gpc.m
%  Przeznaczenie: algorytm w postaci funkcji do celów porównawczych. 
% -------------------------------------------------------------------------

function [y_zad_h, y_zad_pH, yp_h, yp_pH, q1, q3, e] = c_gpc(N, Nu, lambda)
    global Ts
    init_reactor_params;
    %% Parametry
    % Start i stop
    kstart = 3;
    kend = 150;

    % Predykcyjne
    Hp = N;
    Hs = Nu;
%     lambda = 0.001;
    wsplambda = lambda*[1 1];
    wspmi = 1*[1 1];

    % Rzêdy
    na = 2;
    nb = na;

    % Konfiguracja modelu
    ny = 2;
    nu = 2;

    % Warunki pocz¹tkowe
    q10 = 16.6;                             
    q30 = 15.6;                             

    x0 = initial_conditions(q10, q30);      
    Wa0 = x0(1);
    Wb0 = x0(2);
    h0 = x0(3);
    pH0 = x0(4);

    %% Model liniowy
    load model_linear.mat;
    %     load model_linear_nlhw.mat;

    %% Obliczenie macierzy
    [G, K] = calculate_matrices(Hp,Hs,wsplambda,wspmi,ny,nu,na,nb);

    %% Alokacja wektorów
    % Trajektora zadana na ca³ym horyzoncie
    y_zad_h = zeros(kend,1);
    y_zad_pH = zeros(kend,1);
    % Trajektoria zadana w danej iteracji algorytmu
    y_ref_h = zeros(Hp,1);
    y_ref_pH = zeros(Hp,1);
    y_ref = zeros(2*Hp,1);
    % Trajektoria wymuszona
    y_pred_h = zeros(Hp,1);
    y_pred_pH = zeros(Hp,1);
    % Wyjœcie procesu
    yp_h = zeros(kend,1);
    yp_pH = zeros(kend,1);
    yp_h(1:kstart-1) = h0;
    yp_pH(1:kstart-1) = pH0;   
    % Wyjœcie modelu
    ym_h = zeros(kend,1);
    ym_pH = zeros(kend,1);
    % Trajektoria swobodna
    y0_h = zeros(Hp,1);
    y0_pH = zeros(Hp,1);
    y0 = zeros(2*Hp,1);
    % Przyrosty sterowan
    detu = zeros(2*Hs,1);
    % Niemierzealne zak³ócenia i niedok³adnoœci modelu
    d_h = zeros(kend,1);
    d_pH = zeros(kend,1);
    % Sterowania modelu
    u = zeros(2,kend);
    % Sterowania procesu
    q1(1:kstart-1) = q10;
    q3(1:kstart-1) = q30;
    %% Wartoœci zadaane
    y_zad_pH(1:kstart) = pH0;
    y_zad_pH(3:19) = 6;
    y_zad_pH(20:39) = 8;
    y_zad_pH(40:59) = 5;
    y_zad_pH(60:79) = 9;
    y_zad_pH(80:99) = 4;
    y_zad_pH(100:kend) = 10;

    y_zad_h(1:kend) = h0;
    
    %% Ograniczenia
    q1min = 0;
    q1max = 30;
    q3min = 0;
    q3max = 30;
    deltaq1min = q1min - q1max;
    deltaq1max = q1max - q1min;
    deltaq3min = q3min - q3max;
    deltaq3max = q3max - q3min;
    u1min = (q1min - q10)/15;
    u1max = (q1max - q10)/15;
    u2min = (q3min - q30)/15;
    u2max = (q3max - q30)/15;
    deltau1min = u1min - u1max;
    deltau1max = u1max - u1min;
    deltau2min = u2min - u2max;
    deltau2max = u2max - u2min;
    
    %% Algorytm GPC
    u(1,1:kstart-1) = (q1(1:kstart-1)-q10)/15;
    u(2,1:kstart-1) = (q3(1:kstart-1)-q30)/15;
    x0 = x0(:,1:3);
    for k = kstart:kend
        %% Symulacja procesu
        % yp_h yp_pH
        [t,x] = ode45(@(t,x) sim_reactor(t, x, q1(k-1), q3(k-1)), [(k-1)*Ts (k)*Ts],x0);
        x0 = [x(end,1) x(end,2) x(end,3)];
        yp_h(k) = x(end,3);
        yp_pH(k) = calc_pH(x(end,1), x(end,2));

        % Skalowanie
        yp_h_sc(k) = (yp_h(k)-h0)/20;
        yp_pH_sc(k) = (yp_pH(k)-pH0)/5;
        %% Symulacja modelu liniowego
        % y1 h
        for n = 1:nu
            for i = 1:nb
                ym_h(k) = ym_h(k) + b(1,n,i)*u(n,k-i);
            end
        end
        for i = 1:na
            ym_h(k) = ym_h(k) - a(1,i)*ym_h(k-i);
        end        
        % y2 pH
        for n = 1:nu
            for i = 1:nb
                ym_pH(k) = ym_pH(k) + b(2,n,i)*u(n,k-i);
            end
        end
        for i = 1:na
            ym_pH(k) = ym_pH(k) - a(2,i)*ym_pH(k-i);
        end        
        % Niepewnoœci pomiaru i modelu
        d_h(k) = yp_h_sc(k) - ym_h(k);
        d_pH(k) = yp_pH_sc(k) - ym_pH(k);
        % Obliczenie trajektori referencyjnej
        y_ref_h = (y_zad_h(k)-h0)*ones(Hp,1)/20;
        y_ref_pH = (y_zad_pH(k)-pH0)*ones(Hp,1)/5;

        % Obliczenie y0 h
        p = 1;
        y0_h(p) = d_h(k)+b(1,1,1)*u(1,k-1)+b(1,1,2)*u(1,k-1)+b(1,2,1)*u(2,k-1)+b(1,2,2)*u(2,k-1)-a(1,1)*ym_h(k)-a(1,2)*ym_h(k-1);

        p = 2;
        y0_h(p) = d_h(k)+b(1,1,1)*u(1,k-1)+b(1,1,2)*u(1,k-1)+b(1,2,1)*u(2,k-1)+b(1,2,2)*u(2,k-1)-a(1,1)*y0_h(p-1)-a(1,2)*ym_h(k);

        for p = 3:Hp
            y0_h(p) = d_h(k)+b(1,1,1)*u(1,k-1)+b(1,1,2)*u(1,k-1)+b(1,2,1)*u(2,k-1)+b(1,2,2)*u(2,k-1)-a(1,1)*y0_h(p-1)-a(1,2)*y0_h(p-2); 
        end
        % Obliczenie y0 pH
        p = 1;
        y0_pH(p) = d_pH(k)+b(2,1,1)*u(1,k-1)+b(2,1,2)*u(1,k-1)+b(2,2,1)*u(2,k-1)+b(2,2,2)*u(2,k-1)-a(2,1)*ym_pH(k)-a(2,2)*ym_pH(k-1); 

        p = 2;
        y0_pH(p) = d_pH(k)+b(2,1,1)*u(1,k-1)+b(2,1,2)*u(1,k-1)+b(2,2,1)*u(2,k-1)+b(2,2,2)*u(2,k-1)-a(2,1)*y0_pH(p-1)-a(2,2)*ym_pH(k); 

        for p = 3:Hp
            y0_pH(p) = d_pH(k)+b(2,1,1)*u(1,k-1)+b(2,1,2)*u(1,k-1)+b(2,2,1)*u(2,k-1)+b(2,2,2)*u(2,k-1)-a(2,1)*y0_pH(p-1)-a(2,2)*y0_pH(p-2);
        end
        % Z³o¿enie wektora y_ref
        for i = 1:Hp
           y_ref(2*i-1) = y_ref_h(i);
           y_ref(2*i) = y_ref_pH(i);
        end
        % Z³o¿enie wektora y0
        for i = 1:Hp
           y0(2*i-1) = y0_h(i);
           y0(2*i) = y0_pH(i);
        end       
        % Obliczenie detu
        detu = K*(y_ref - y0);
        % Obliczenie u
        if detu(1) > deltau1max
            detu(1) = deltau1max;
        end
        if detu(1) < deltau1min
            detu(1) = deltau1min;
        end
        if detu(2) > deltau2max
            detu(2) = deltau2max;
        end
        if detu(2) < deltau2min
            detu(2) = deltau2min;
        end
        u(1,k) = detu(1) + u(1,k-1);
        u(2,k) = detu(2) + u(2,k-1);
        if u(1,k) > u1max
            u(1,k) = u1max;
        end
        if u(1,k) < u1min
            u(1,k) = u1min;
        end
        if u(2,k) > u2max
            u(2,k) = u2max;
        end
        if u(2,k) < u2min
            u(2,k) = u2min;
        end  
        q1(k) = u(1,k)*15 + q10;
        q3(k) = u(2,k)*15 + q30;      
    end
    e(1)=(y_zad_h - yp_h)'*(y_zad_h - yp_h);
    e(2)=(y_zad_pH - yp_pH)'*(y_zad_pH - yp_pH);
end
