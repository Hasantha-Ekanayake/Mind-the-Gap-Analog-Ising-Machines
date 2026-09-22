clc;
clear all;
close all;


nOsc=15;
J = load('random_15.csv');

% ------------------Model-------------


K=1;
Ac = K;
An = 0.001;
tstop = 40; tstep = 2e-3;
 
 a1.k =  4/ tstop; % Compute k
 f1 = @(t, args)  t * args.k;

 F1 = @(t,X) KuramotoF1_DIM(X, Ac, f1(t,a1), nOsc,J);
 G1 = @(t,X) An*eye(nOsc);
 
 X = 0.5*ones(nOsc,1); 
 obj1 = sde(F1, G1, 'StartState', X);
 [S1, T1] = simulate(obj1, tstop/tstep, 'DeltaTime', tstep);
 
    
 figure(1); 
 
 % Evaluate f1(t) over the time vector
 f1_values = arrayfun(@(t_val) f1(t_val, a1), T1); % Apply f1 to each time value

 % Add a constant K = 1 line
 yline(1, '-', 'Color', "#0072BD", 'LineWidth', 4);
 hold on;
 plot(T1, f1_values, 'LineWidth', 4, 'Color', "#D95319"); 

 figure(2); 
 plot(f1_values,S1, LineWidth=4);
 hold on
 xline(30, '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]); % Gray line
 hold on
 
 xlabel('K_S');xlim([0, max(f1_values)]);
 ylabel('phase(\pi)');

 for k = 1:length(T1)
  ix = find(mod(round(S1(k,:)), 2));
  cuts1(k) = -sum(sum(J(ix, setdiff(1:nOsc, ix))));
 end 
 MC =max(cuts1)
 % Adjust axes
 ax = gca;
 ax.LineWidth = 3; % Thicker axes
 ax.FontSize = 20;   % Larger tick labels 

% -------------DIM functions----------------------------


function fout = KuramotoF1_DIM(x, Ac, As, n, J)

for c = 1:n

      fout(c, 1) = - Ac * J(c, :) * sin(pi*(x(c) + x));
      
end

fout = (fout-As*sin(2*pi*x))/pi;

end

