%% ------------------------------------------------------------------------
%                       Praca dyplomowa magisterska
%  Autor: Patryk Okoñski
%  Temat: Nieliniowa regulacja predykcyjna procesu dwuwymiarowego
%  wykorzystuj¹ca neuronowy model Wienera.
%  ------------------------------------------------------------------------

%% ---------------------- FILE INFO ---------------------------------------
%  Nazwa: c_mpc_npsl.m
%  Przeznaczenie: algorytm w postaci funkcji do celów porównawczych. 
% -------------------------------------------------------------------------

function [y_zad_h, y_zad_pH, yp_h, yp_pH, q1, q3, e] = c_mpc_npsl(N, Nu, lambda)

    global Ts

    init_reactor_params;

    %% Neuronowe modele wienera
    load('Dane/model_h.mat');
    load('Dane/model_pH.mat');

    % Rzêdy
    nb = max(model_h.nb)+1;
    nf = max(model_h.nf)+1;
    nk = max(model_h.nk)+1;
    % Konfiguracja modelu
    ny = 2;
    nu = 2;    
    na = nb;
    %% Parametry
    % Regulacji
    kstart = 3;
    kend = 150;

    % Predykcyjne
    Hp = N;
    Hs = Nu;
%     lambda = 0.1;
    wsplambda = lambda*[1 1];
    wspmi = 1*[1 1];

    % Warunki pocz¹tkowe
    q10 = 16.6;                             
    q30 = 15.6;                             

    x0 = initial_conditions(q10, q30);      
    Wa0 = x0(1);
    Wb0 = x0(2);
    h0 = x0(3);
    pH0 = x0(4);
    x0 = x0(:,1:3);

    % -----Parametry neuronowego modelu wienera na potrzeby symulacji
    % Model wyjœcia h
    rnn_h = model_h.OutputNonlinearity.Parameters.RegressorMean;
    Qnn_h = model_h.OutputNonlinearity.Parameters.NonLinearSubspace;
    Pnn_h = model_h.OutputNonlinearity.Parameters.LinearSubspace;
    Lnn_h = model_h.OutputNonlinearity.Parameters.LinearCoef;
    lnn_h = Pnn_h*Lnn_h;
    bnn_h = model_h.OutputNonlinearity.Parameters.Dilation;
    cnn_h = model_h.OutputNonlinearity.Parameters.Translation;
    ann_h = model_h.OutputNonlinearity.Parameters.OutputCoef;
    dnn_h = model_h.OutputNonlinearity.Parameters.OutputOffset;
    bnn_h = Qnn_h*bnn_h;

    B11_h = model_h.B{1}(2); B21_h = model_h.B{2}(2);
    F11_h = model_h.F{1}(2); F21_h = model_h.F{2}(2);
    aw_h(1) = (F11_h + F21_h); aw_h(2) = F11_h*F21_h;
    bw_h(1,1) = B11_h; bw_h(1,2) = B11_h*F21_h;
    bw_h(2,1) = B21_h; bw_h(2,2) = B21_h*F11_h;

    % Model wyjœcia pH
    rnn_pH = model_pH.OutputNonlinearity.Parameters.RegressorMean;
    Qnn_pH = model_pH.OutputNonlinearity.Parameters.NonLinearSubspace;
    Pnn_pH = model_pH.OutputNonlinearity.Parameters.LinearSubspace;
    Lnn_pH = model_pH.OutputNonlinearity.Parameters.LinearCoef;
    lnn_pH = Pnn_pH*Lnn_pH;
    bnn_pH = model_pH.OutputNonlinearity.Parameters.Dilation;
    cnn_pH = model_pH.OutputNonlinearity.Parameters.Translation;
    ann_pH = model_pH.OutputNonlinearity.Parameters.OutputCoef;
    dnn_pH = model_pH.OutputNonlinearity.Parameters.OutputOffset;
    bnn_pH = Qnn_pH*bnn_pH;

    B11_pH = model_pH.B{1}(2); B21_pH = model_pH.B{2}(2);
    F11_pH = model_pH.F{1}(2); F21_pH = model_pH.F{2}(2);
    aw_pH(1) = (F11_pH + F21_pH); aw_pH(2) = F11_pH*F21_pH;
    bw_pH(1,1) = B11_pH; bw_pH(1,2) = B11_pH*F21_pH;
    bw_pH(2,1) = B21_pH; bw_pH(2,2) = B21_pH*F11_pH;	

    Kn = model_h.OutputNonlinearity.NumberOfUnits;
    % Zmiana formatu parametrów a i b
    b(1,1,1) = bw_h(1,1);
    b(1,1,2) = bw_h(1,2);
    b(1,2,1) = bw_h(2,1);
    b(1,2,2) = bw_h(2,2);
    b(2,1,1) = bw_pH(1,1);
    b(2,1,2) = bw_pH(1,2);
    b(2,2,1) = bw_pH(2,1);
    b(2,2,2) = bw_pH(2,2);
    a(1,1) = aw_h(1);
    a(1,2) = aw_h(2);
    a(2,1) = aw_pH(1);
    a(2,2) = aw_pH(2);

    %% Alokacja wektorów
    % Trajektora zadana na ca³ym horyzoncie symulacji
    y_zad_h = zeros(kend,1);
    y_zad_pH = zeros(kend,1);
    % Trajektoria zadana w danej iteracji algorytmu
    y_ref_h = zeros(Hp,1);
    y_ref_pH = zeros(Hp,1);
    y_ref = zeros(2*Hp,1);
    % Wyjœcia procesu
    yp_h = zeros(kend,1);
    yp_pH = zeros(kend,1);
    yp_h(1:kstart-1) = h0;
    yp_pH(1:kstart-1) = pH0;
    % Skalowane wyjœcia procesu
    yp_h_sc(1:kstart-1) = (yp_h(1:kstart-1)-h0)/20;
    yp_pH_sc(1:kstart-1) = (yp_pH(1:kstart-1)-pH0)/5;
    % Wyjœcia modelu
    ym_h = zeros(kend,1);
    ym_pH = zeros(kend,1);
    % Trajektoria swobodna
    y0_h = zeros(Hp,1);
    y0_pH = zeros(Hp,1);
    y0 = zeros(2*Hp,1);
    % Przyrosty sterowañ
    detu = zeros(2*Hs,1);
    % Niemierzealne zak³ócenia i niedok³adnoœci
    d_h = zeros(kend,1);
    d_pH = zeros(kend,1);
    % Sterowania procesu
    q1(1:kstart-1) = q10;
    q3(1:kstart-1) = q30;
    % Sterowania modelu
    u = zeros(2,kend);
    u(1,1:kstart-1) = (q1(1:kstart-1)-q10)/15;
    u(2,1:kstart-1) = (q3(1:kstart-1)-q30)/15;
    %% Wartoœci zadane
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
    
    % Dodatkowe sygna³y na potrzeby symulacji modelu
    v_h(1:kstart-1) = 0;
    v_pH(1:kstart-1) = 0;
    v_h0(1:Hp) = 0;
    v_pH0(1:Hp) = 0;

    %% ALGORYTM MPC-NPSL
    % Obliczenie macierzy S. Format S(:,:,j) = [s(1,1,j) s(1,2,j); s(2,1,j) s(2,2,j)]   s(m,n,j) = wartosc;
    s = zeros(2,2,Hp);
    S = zeros(2,2);
    Sp = [];
    G = [];
    K = [];
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
    % Macierze MI i LAMBDA na potrzeby wyznaczenia macierzy K
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
    % Pêtla g³ówna
    for k = kstart:kend
        % -------------------------------------- Symulacja procesu
        % yp_h yp_pH
        [t,x] = ode45(@(t,x) sim_reactor(t, x, q1(k-1), q3(k-1)), [(k-1)*Ts (k)*Ts],x0);
        x0 = [x(end,1) x(end,2) x(end,3)];
        yp_h(k) = x(end,3);
        yp_pH(k) = calc_pH(x(end,1), x(end,2));

        % Scaling
        yp_h_sc(k) = (yp_h(k)-h0)/20;
        yp_pH_sc(k) = (yp_pH(k)-pH0)/5;

        % -------------------------------------- Symulacja modelu
        % ---------- Obliczenie wyjœcia h
        v_h(k) = bw_h(1,1)*u(1,k-1) + bw_h(1,2)*u(1,k-2)...
        + bw_h(2,1)*u(2,k-1) + bw_h(2,2)*u(2,k-2)...
        - aw_h(1)*v_h(k-1) - aw_h(2)*v_h(k-2);

        ym_h(k) = lnn_h*(v_h(k) - rnn_h) + dnn_h;

        for i = 1:Kn
            ym_h(k) = ym_h(k) + ann_h(i)*logsig((v_h(k) - rnn_h)*bnn_h(i) + cnn_h(i));
        end
        % ---------- Obliczenie wyjœcia pH
        v_pH(k) = bw_pH(1,1)*u(1,k-1) + bw_pH(1,2)*u(1,k-2)...
        + bw_pH(2,1)*u(2,k-1) + bw_pH(2,2)*u(2,k-2)...
        - aw_pH(1)*v_pH(k-1) - aw_pH(2)*v_pH(k-2);

        ym_pH(k) = lnn_pH*(v_pH(k) - rnn_pH) + dnn_pH;

        for i = 1:Kn
            ym_pH(k) = ym_pH(k) + ann_pH(i)*logsig((v_pH(k) - rnn_pH)*bnn_pH(i) + cnn_pH(i));
        end
        % -------------------------------------- Inne obliczenia    
        % Niedok³adnoœci i niepewnoœci pomiaru
        d_h(k) = yp_h_sc(k) - ym_h(k);
        d_pH(k) = yp_pH_sc(k) - ym_pH(k);
        % Trajektoria referencyjna
        y_ref_h = (y_zad_h(k)-h0)*ones(Hp,1)/20;
        y_ref_pH = (y_zad_pH(k)-pH0)*ones(Hp,1)/5;

        % -------------------------------------- Obliczenia dla MPC-NPSL
        % Obliczenie Kw
        %-----k11, k12
        k11 = 0;
        for i = 1:Kn
            zh = (v_h(k) - rnn_h)*bnn_h(i) + cnn_h(i);
            k11 = k11 + ann_h(i)*logsig(zh)*(1-logsig(zh))*bnn_h(i);
        end
        k12 = k11;
        %-----k21, k22
        k21 = 0;
        for i = 1:Kn
            zpH = (v_pH(k) - rnn_pH)*bnn_pH(i) + cnn_pH(i);
            k21 = k21 + ann_pH(i)*logsig(zpH)*(1-logsig(zpH))*bnn_pH(i);
        end
        k22 = k21;

        % Drbugowanie k
        k11k(k) = k11;
        k21k(k) = k21;

        Kw = [k11 k12; k21 k22];

        % Obliczenie Sp
        for j = 1:Hp
            Sp(:,:,j) = Kw.*S(:,:,j);
        end


        % Obliczenie G
        for i = 1:Hs
            for j = i:Hp
                G((2*j-1):2*j,(2*i-1):2*i) = Sp(:,:,j-i+1);
            end
        end

        % Obliczenie K
        K=inv(G'*MI*G+LAMBDA)*G'*MI;

        % Obliczenie v0
        p = 1;
        v_h0(p) = b(1,1,1)*u(1,k-1)+b(1,1,2)*u(1,k-1)+b(1,2,1)*u(2,k-1)+b(1,2,2)*u(2,k-1)-a(1,1)*v_h(k)-a(1,2)*v_h(k-1);
        v_pH0(p) = b(2,1,1)*u(1,k-1)+b(2,1,2)*u(1,k-1)+b(2,2,1)*u(2,k-1)+b(2,2,2)*u(2,k-1)-a(2,1)*v_pH(k)-a(2,2)*v_pH(k-1);

        p = 2;
        v_h0(p) = b(1,1,1)*u(1,k-1)+b(1,1,2)*u(1,k-1)+b(1,2,1)*u(2,k-1)+b(1,2,2)*u(2,k-1)-a(1,1)*v_h0(1)-a(1,2)*v_h(k);		
        v_pH0(p) = b(2,1,1)*u(1,k-1)+b(2,1,2)*u(1,k-1)+b(2,2,1)*u(2,k-1)+b(2,2,2)*u(2,k-1)-a(2,1)*v_pH0(1)-a(2,2)*v_pH(k);

        for p = 3:Hp
            v_h0(p) = b(1,1,1)*u(1,k-1)+b(1,1,2)*u(1,k-1)+b(1,2,1)*u(2,k-1)+b(1,2,2)*u(2,k-1)-a(1,1)*v_h0(p-1)-a(1,2)*v_h0(p-2);
            v_pH0(p) = b(2,1,1)*u(1,k-1)+b(2,1,2)*u(1,k-1)+b(2,2,1)*u(2,k-1)+b(2,2,2)*u(2,k-1)-a(2,1)*v_pH0(p-1)-a(2,2)*v_pH0(p-2);
        end

        % Obliczenie y0
        for p = 1:Hp
            y0_h(p) = dnn_h + lnn_h*(v_h0(p) - rnn_h) + d_h(k);
            y0_pH(p) = dnn_pH + lnn_pH*(v_pH0(p) - rnn_pH) + d_pH(k);
            for i = 1:Kn
                y0_h(p) = y0_h(p) + ann_h(i)*logsig(bnn_h(i)*(v_h0(p) - rnn_h) + cnn_h(i));
                y0_pH(p) = y0_pH(p) + ann_pH(i)*logsig(bnn_pH(i)*(v_pH0(p) - rnn_pH) + cnn_pH(i));
            end
        end

        % Budowanie wektora y_ref
        for i = 1:Hp
           y_ref(2*i-1) = y_ref_h(i);
           y_ref(2*i) = y_ref_pH(i);
        end

        % Budowanie wektora y0
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

%     figure(1)
%     grid on;
%     subplot(2,1,1);
%     stairs(q1);
%     ylabel('q_1');
%     xlabel('Krok [k]');
%     subplot(2,1,2);
%     stairs(q3);
%     ylabel('q_3');
%     xlabel('Krok [k]');
% 
%     figure(2)
%     grid on;
%     subplot(2,1,1);
%     stairs(y_zad_h,'r');
%     hold on;
%     plot(yp_h,'--b');
%     hold off
%     ylabel('h_{zad}, h');
%     xlabel('Krok [k]');
%     legend('h_{zad}','h')
%     title(['Regulacja MPC-NPSL [h], MSE_h=',num2str(e1)]);
% 
%     subplot(2,1,2);
%     stairs(y_zad_pH,'r');
%     hold on;
%     plot(yp_pH,'--b');
%     hold off
%     ylabel('pH_{zad}, pH');
%     xlabel('Krok [k]');
%     legend('pH_{zad}','pH');
%     title(['Regulacja MPC-NPSL [pH], MSE_{pH}=',num2str(e2)]);

    % figure(3)
    % grid on
    % stairs(k11k)
    % hold on
    % stairs(k21k,'--r')
    % hold off
    % xlabel('k')
    % ylabel('k11, k21')
    % legend('k11k','k21k')
end