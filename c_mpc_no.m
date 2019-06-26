%% ------------------------------------------------------------------------
%                       Praca dyplomowa magisterska
%  Autor: Patryk Okoñski
%  Temat: Nieliniowa regulacja predykcyjna procesu dwuwymiarowego
%  wykorzystuj¹ca neuronowy model Wienera.
%  ------------------------------------------------------------------------

%% ---------------------- FILE INFO ---------------------------------------
%  Nazwa: c_mpc_no.m
%  Przeznaczenie: algorytm w postaci funkcji do celów porównawczych. 
% -------------------------------------------------------------------------

function [yref_h, yref_pH, ymod_h, ymod_pH, q1, q3, e] = c_mpc_no(Na, Nua, lambdaa)
    
    global Ts
    global N Nu lambda k u u1 u2
    global v_h bw_h aw_h lnn_h rnn_h dnn_h d_h ann_h cnn_h bnn_h K
    global v_pH bw_pH aw_pH lnn_pH rnn_pH dnn_pH d_pH ann_pH cnn_pH bnn_pH
    global ymod_h_sim ymod_pH_sim yref_h_sim yref_pH_sim
    
    N = Na;
    Nu = Nua;
    lambda = lambdaa;
    
    init_reactor_params;

    %% Neuronowe modele wienera wyjœæ h i pH
    load('Dane/model_h.mat');
    load('Dane/model_pH.mat');

    % Rzêdy
    nb = max(model_h.nb);
    nf = max(model_h.nf);
    nk = max(model_h.nk);

    %% USTAWIENIA
    % MPC-NO 
%     N = 10;
%     Nu = 3;
%     lambda = 0.1;

    q1min = 0;
    q1max = 30;
    q3min = 0;
    q3max = 30;
    deltaq1min = q1min - q1max;
    deltaq1max = q1max - q1min;
    deltaq3min = q3min - q3max;
    deltaq3max = q3max - q3min;

    % Warunki pocz¹tkowe
    q10 = 16.6;                             
    q30 = 15.6;                             

    x0 = initial_conditions(q10, q30);      
    Wa0 = x0(1);
    Wb0 = x0(2);
    h0 = x0(3);
    pH0 = x0(4);    

    u1min = (q1min - q10)/15;
    u1max = (q1max - q10)/15;
    u2min = (q3min - q30)/15;
    u2max = (q3max - q30)/15;    

    umin = [u1min u2min];
    umax = [u1max u2max];

    deltau1min = u1min - u1max;
    deltau1max = u1max - u1min;
    deltau2min = u2min - u2max;
    deltau2max = u2max - u2min;     

    deltaumin = [deltau1min deltau2min];
    deltaumax = [deltau1max deltau2max];

    kmax = 150;
    kmin = 3;

    %% Alokacja wektorów
    % -----Wartoœci zadane
    % h   
    yref_h(1:kmax-N) = h0;
    % pH
    yref_pH(1:2) = pH0;
    yref_pH(3:19) = 6;
    yref_pH(20:39) = 8;
    yref_pH(40:59) = 5;
    yref_pH(60:79) = 9;
    yref_pH(80:99) = 4;
    yref_pH(100:kmax-N) = 10;

    % Sterowania procesu
    q1(1:kmin-1) = q10;
    q3(1:kmin-1) = q30;
    % Wyjœcia procesu
    ymod_h(1:kmin-1) = h0;
    ymod_pH(1:kmin-1) = pH0;

    % Sterowania i wyjœcia modelu
    u1 = zeros(1,kmax+Nu-1);
    u2 = zeros(1,kmax+Nu-1);
    ymod_h_sim = zeros(1,kmax+N);
    ymod_pH_sim = zeros(1,kmax+N);

    u1(1:kmin-1) = (q1(1:kmin-1)-q10)/15;
    u2(1:kmin-1) = (q3(1:kmin-1)-q30)/15;
    u = [u1; u2];
    ymod_h_sim(1:kmin-1) = (ymod_h(1:kmin-1)-h0)/20;
    ymod_pH_sim(1:kmin-1) = (ymod_pH(1:kmin-1)-pH0)/5;

    % Dodatkowe sygna³y na potrzeby symulacji modelu
    v_h(1:kmin-1) = 0;
    v_pH(1:kmin-1) = 0;  

    % Trajektoria referencyjna
    yref_pH_sim = (yref_pH-pH0)/5;
    yref_h_sim = (yref_h-h0)/20;

    yref_sp_vec = zeros(2*N,1);

    % -----Parametry na potrzeby symulacji neuronowego modelu wienera
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

    K = model_h.OutputNonlinearity.NumberOfUnits;

    options_fmincon = optimoptions(@fmincon,'Algorithm','sqp', ...
    'Display','none');

    du_opt_init = zeros(2*Nu,1);

    A_opt = [-kron(tril(ones(Nu,Nu),0),eye(2)); ...
        kron(tril(ones(Nu,Nu),0),eye(2))];
    lb_opt = kron(ones(1,Nu),deltaumin);
    ub_opt = kron(ones(1,Nu),deltaumax);    

    x0 = x0(:,1:3);
    tic
    for k = kmin:kmax-N
        % -------------------------------------- Symulacja procesu
        [t,x] = ode45(@(t,x) sim_reactor(t, x, q1(k-1), q3(k-1)), [(k-1)*Ts (k)*Ts],x0);
        x0 = [x(end,1) x(end,2) x(end,3)];
        ymod_h(k) = x(end,3);
        ymod_pH(k) = calc_pH(x(end,1), x(end,2));

        % Skalowanie
        ymod_h_sim(k) = (ymod_h(k)-h0)/20;
        ymod_pH_sim(k) = (ymod_pH(k)-pH0)/5;


        % -------------------------------------- Symulacja modelu
        % ---------- Obliczenie wyjœcia h
        v_h(k) = bw_h(1,1)*u1(k-1) + bw_h(1,2)*u1(k-2)...
            + bw_h(2,1)*u2(k-1) + bw_h(2,2)*u2(k-2)...
            - aw_h(1)*v_h(k-1) - aw_h(2)*v_h(k-2);

        yh = lnn_h*(v_h(k) - rnn_h) + dnn_h;

        for i = 1:K
            yh = yh + ann_h(i)*logsig((v_h(k) - rnn_h)*bnn_h(i) + cnn_h(i));
        end
        % ---------- Obliczenie wyjœcia pH
        v_pH(k) = bw_pH(1,1)*u1(k-1) + bw_pH(1,2)*u1(k-2)...
            + bw_pH(2,1)*u2(k-1) + bw_pH(2,2)*u2(k-2)...
            - aw_pH(1)*v_pH(k-1) - aw_pH(2)*v_pH(k-2);

        ypH = lnn_pH*(v_pH(k) - rnn_pH) + dnn_pH;

        for i = 1:K
            ypH = ypH + ann_pH(i)*logsig((v_pH(k) - rnn_pH)*bnn_pH(i) + cnn_pH(i));
        end

        % Niepewnoœci modelu i niemierzalne zak³ócenia
        d_h(k) = ymod_h_sim(k) - yh;
        d_pH(k) = ymod_pH_sim(k) - ypH;

        % 
        b_opt = [-kron(ones(1,Nu),umin)+kron(ones(1,Nu),u(:,k-1)') ...
        kron(ones(1,Nu),umax)-kron(ones(1,Nu),u(:,k-1)')];        


        du_opt = fmincon(@model_mpc, du_opt_init, A_opt, b_opt, [], [], ...
            lb_opt, ub_opt, [], options_fmincon);

        u(1:2,k) = du_opt(1:2)+u(1:2,k-1);
        u1(k) = u(1,k);
        u2(k) = u(2,k);

        % Skalowanie sterowañ
        q1(k) = u1(k)*15 + q10;
        q3(k) = u2(k)*15 + q30;
    end
    czas = toc;
    ct = max(size(ymod_h));
    e(1) = (yref_h(1:ct) - ymod_h)*(yref_h(1:ct) - ymod_h)';
    e(2) = (yref_pH(1:ct) - ymod_pH)*(yref_pH(1:ct) - ymod_pH)';

    %% Dane do strojenie N
    % save (['Dane/Strojenie/mpcno_strojenie_n_',num2str(N),'.mat'],'yref_h','yref_pH','ymod_h','ymod_pH','q1','q3','e','czas');
    %% Dane do strojenia Nu
    % save (['Dane/Strojenie/mpcno_strojenie_nu_',num2str(Nu),'.mat'],'yref_h','yref_pH','ymod_h','ymod_pH','q1','q3','e','czas');
    %% Dane do strojenia lambda
    % save (['Dane/Strojenie/mpcno_strojenie_lambda_',num2str(lambda),'.mat'],'yref_h','yref_pH','ymod_h','ymod_pH','q1','q3','e','czas');


%     w1 = figure(1)
%     subplot(2,1,1)
%     grid on;
%     stairs(yref_h,'r');
%     hold on
%     plot(ymod_h,'--');  
%     hold off
%     xlabel('Krok [k]');
%     ylabel('h [m]');
%     legend('h_{zad}','h');
%     title(['Regulacja MPC-NO [h], MSE_{h}=',num2str(e(1))]);
% 
%     subplot(2,1,2)
%     grid on;
%     stairs(yref_pH,'r');
%     hold on
%     plot(ymod_pH,'--');
%     hold off
%     xlabel('Krok [k]');
%     ylabel('pH [-]');
%     legend('pH_{zad}','pH');
%     title(['Regulacja MPC-NO [pH], MSE_{pH}=',num2str(e(2))]);
% 
%     w2 = figure(2);
%     subplot(2,1,1)
%     grid on;
%     stairs(q1);
%     xlabel('Krok [k]');
%     ylabel('q_1');
%     title('Sterowanie q_1');
% 
%     subplot(2,1,2)
%     grid on;
%     stairs(q3);
%     xlabel('Krok [k]');
%     ylabel('q_3');
%     title('Sterowanie q_3');

    % save ('Dane/mpcno_data.mat');
    % saveas(w1, 'Dane/Wykresy/mpc_no_h','png');
    % saveas(w2, 'Dane/Wykresy/mpc_no_pH','png');
end

