%% ------------------------------------------------------------------------
%                       Praca dyplomowa magisterska
%  Autor: Patryk Okoñski
%  Temat: Nieliniowa regulacja predykcyjna procesu dwuwymiarowego
%  wykorzystuj¹ca neuronowy model Wienera.
%  ------------------------------------------------------------------------

%% ---------------------- FILE INFO ---------------------------------------
%  Nazwa: model_mpc.m
%  Przeznaczenie: Funkcja kosztu na potrzeby MPC-NO.
% -------------------------------------------------------------------------

function [J] = model_mpc(x)
    global Ts
    global N Nu lambda k u u1 u2
    global v_h bw_h aw_h lnn_h rnn_h dnn_h d_h ann_h cnn_h bnn_h K
    global v_pH bw_pH aw_pH lnn_pH rnn_pH dnn_pH d_pH ann_pH cnn_pH bnn_pH
    global ymod_h_sim ymod_pH_sim yref_h_sim yref_pH_sim

    vector_u = kron(tril(ones(Nu,Nu),0),eye(2))*x ...
        +kron(ones(Nu,1),u(:,k-1)); 
    
    u(1,k:k+Nu-1) = vector_u(1:2:end);
    u(2,k:k+Nu-1) = vector_u(2:2:end);
    u(1,k+Nu:k+N) = u(1,k+Nu-1);
    u(2,k+Nu:k+N) = u(2,k+Nu-1);
    
    for p = 1:N 
        % ---------- obliczenie wyjœcia H
		v_h(k+p) = bw_h(1,1)*u(1,k+p-1) + bw_h(1,2)*u(1,k+p-2)...
			+ bw_h(2,1)*u(2,k+p-1) + bw_h(2,2)*u(2,k+p-2)...
			- aw_h(1)*v_h(k+p-1) - aw_h(2)*v_h(k+p-2);
		
		ymod_h_sim(k+p) = lnn_h*(v_h(k+p) - rnn_h) + dnn_h + d_h(k);
		for i = 1:K
			ymod_h_sim(k+p) = ymod_h_sim(k+p) + ann_h(i)*logsig((v_h(k+p) - rnn_h)*bnn_h(i) + cnn_h(i));
        end
        
		% ---------- obliczenie wyjœcia pH
		v_pH(k+p) = bw_pH(1,1)*u(1,k+p-1) + bw_pH(1,2)*u(1,k+p-2)...
			+ bw_pH(2,1)*u(2,k+p-1) + bw_pH(2,2)*u(2,k+p-2)...
			- aw_pH(1)*v_pH(k+p-1) - aw_pH(2)*v_pH(k+p-2);
		
		ymod_pH_sim(k+p) = lnn_pH*(v_pH(k+p) - rnn_pH) + dnn_pH + d_pH(k);

		for i = 1:K
			ymod_pH_sim(k+p) = ymod_pH_sim(k+p) + ann_pH(i)*logsig((v_pH(k+p) - rnn_pH)*bnn_pH(i) + cnn_pH(i));
        end
		
    end
    
   % Funkcja kosztu
    J = norm(yref_h_sim(k)*ones(1,N) - ymod_h_sim(k+1:k+N))^2 ...
        +norm(yref_pH_sim(k)*ones(1,N) - ymod_pH_sim(k+1:k+N))^2 ...
        +lambda*norm(x(1:2:end))^2 ...
        +lambda*norm(x(2:2:end))^2;

end