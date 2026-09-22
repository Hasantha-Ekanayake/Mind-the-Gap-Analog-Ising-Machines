%% 
clc;
clear all;
close all;

nOsc=15;


J = load("random_15.csv");

% ------------------Model-------------




eps = 0.05;    
delta = 1;
Ke = 1;
An = 0.001;
tstop = 40; tstep = 2e-3;
    

 a1.k =  4/ tstop; % Compute k
 f1 = @(t, args)  t * args.k;
 X = zeros(nOsc,1);


 F1 = @(t,X) dynamics_to_sbm(X, f1(t,a1),eps, nOsc,J, delta, Ke);
 G1 = @(t,X) An*eye(nOsc);

 obj1 = sde(F1, G1, 'StartState', X);
 [S1, T1] = simulate(obj1, tstop/tstep, 'DeltaTime', tstep);
     
 for k = 1:length(T1)
        spins      = sign(S1(k,:));   
        ix         = find(spins == -1);
        cuts1(k)   = -sum(sum(J(ix, setdiff(1:nOsc, ix))));

 end

 MC =max(cuts1)



% --------------------Figures---------------------------
    
% --------------------Figures---------------------------

figure(1); 
    

    % Evaluate f1(t) over the time vector
    f1_values = arrayfun(@(t_val) f1(t_val, a1), T1); % Apply f1 to each time value

    % Add a constant K = 1 line
    yline(1, '-', 'Color', "#0072BD", 'LineWidth', 4);
    hold on;
    plot(T1, f1_values, 'LineWidth', 4, 'Color', "#D95319");


    
figure(2)
    plot(f1_values,S1, LineWidth=4);

    ax = gca;
    ax.LineWidth = 3; % Thicker axes
    ax.FontSize = 20;   % Larger tick labels

    xlabel('P'); xlim([0, max(f1_values)]);
    ylabel('x');
    xline(30, '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]); % Gray line



% -------------SBM functions----------------------------


function dxdt = dynamics_to_sbm(x, P, eps, n, J, delta, Ke)

%   x   — n×1 vector of amplitudes x_i
%   K   — coupling strength
%   Ks  — 2nd‐harmonic injection strength
%   eps — scaling parameter ε
%   n   — number of nodes
%   J   — coupling matrix 

  dxdt = zeros(n,1);

  for i = 1:n
    % sum of J over j ≠ i
    sumJ = sum(J(i,:));      % assumes J(i,i)=0
    
    % coupling term: sum_j J_ij * x_j
    sumJx = J(i,:) * x;
    
    % linear growth term
    lin = (-delta+P) * x(i);
    
    % cubic saturation
    cub = - Ke*x(i)^3;
    
    % cross‐injection from other nodes
    coup = eps * sumJx;
    
    % total derivative
    dxdt(i) = lin + cub + coup;
  end
end

