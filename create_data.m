%% 1. Detector noise
% Initial amounts set to 0 to generate only detector noise

% Run simulation
[times1, detector1, logdet1] = simulate_radiation(uint64(0), uint32(0));

% Save as matlab data
save("data1.mat", "times1", "detector1", "logdet1");

% Save as csv
writematrix([times1, detector1, logdet1], 'data1.csv')


%% 2. Secular Equilibrium
% Initial amounts set to values predicted by secular equilibrium

% Calculate lambda values
tp = 30.17*525948.766*60; % T1/2 for 137Cs (in s)
lambda_p = log(2)/tp; % Decay rate constant for 137Cs (s^-1)

td = 2.552*60; % T1/2 for 137mBa (in s)
lambda_d = log(2)/td; % Decay rate constant for 137mBa (s^-1)


% Run simulation
[times2, detector2, logdet2] = simulate_radiation(uint64(5e14), uint32(5e14 * lambda_p/lambda_d));

% Save as matlab data
save("data2.mat", "times2", "detector2", "logdet2");

% Save as csv
writematrix([times2, detector2, logdet2], 'data2.csv')


%% 3. 137mBa Decay
% No 137Cs, only 137mBa.  Simulates decay of 137mBa isolated from parent
% isotope.

% Run simulation
[times3, detector3, logdet3] = simulate_radiation(uint64(0), uint32(5e14*lambda_p/lambda_d));

% Save as matlab data
save("data3.mat", "times3", "detector3", "logdet3");

% Save as csv
writematrix([times3, detector3, logdet3], 'data3.csv')



%% 4. Buildup of 137mBa
% Initial sample contains mostly 137Cs, with only a small amount of 137mBa.
% Simulates decay of a fresh sample of 137Cs.

% Run simulation
[times4, detector4, logdet4] = simulate_radiation(uint64(5e14), uint32(1e7));

% Save as matlab data
save("data4.mat", "times4", "detector4", "logdet4");

% Save as csv
writematrix([times4, detector4, logdet4], 'data4.csv')


