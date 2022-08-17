
dt = 1; % Time step = 1 s; you can change this if required
nt = 2000; % Keep simulation time dt*nt = 2000 s

tp = 30.17*525948.766*60; % T1/2 for 137Cs (in s)
lambda_p = log(2)/tp; % Decay rate constant for 137Cs (s^-1)

td = 2.552*60; % T1/2 for 137mBa (in s)
lambda_d = log(2)/td; % Decay rate constant for 137mBa (s^-1)

n0cs = uint64(5e14); % Initial number of 137Cs nuclei
ncs = n0cs; % Current number of 137Cs nuclei
p1 = 1 - exp(-lambda_p*dt); % Probability of decay of a given 137Cs nucleus per time step

n0ba = uint32(1e7); % Initial number of 137mBa nuclei
nba = n0ba; % Current number of 137mBa nuclei
p2 = 1 - exp(-lambda_d*dt); % Probability of decay of a given 137mBa nucleus per time step

eff = 0.02; % Detector efficiency (2%)
pnoise = 0.8; % Empirical noise parameters
noise = 10/pnoise;

times = [dt:dt:nt*dt]; % Array with the time axis values
detector = 0*times; % Array with detector reads for each time step
logdet = detector; % Log (ln) of the detector reads

for ii = 1:nt
    % The number of 137Cs decays in the current time step.
    % Matlab's binomial RN generator can't handle very large integers,
    % so we approximate a binomial variate by using Gaussian noise.
    % Since the number of decays must be non-negative,
    % any negative values are replaced with 0.
    decays1 = double(ncs)*p1 + normrnd(0, sqrt(double(ncs)*p1*(1-p1)));
    decays1 = max(decays1, 0);

    % Each 137Cs decay decreases the number of 137Cs nuclei
    % and increases the number of 137mBa nuclei:
    ncs = ncs - decays1;
    nba = nba + decays1;

    % The number of 137mBa decays in the current time step.
    % For 137mBa decay we use a true binomial variate when nba is small,
    % otherwise use a Gaussian approximation (which is much faster):
    if (nba < 10000)
        decays2 = binornd(double(nba),p2);
    else
        decays2 = double(nba)*p2 + normrnd(0, sqrt(double(nba)*p2*(1-p2)));
        decays2 = max(decays2, 0);
    end % End of If

    % Each 137mBa decay decreases the number of 137mBa nuclei:
    nba = nba - decays2;
    ii; % Lazy progress indicator
    % Detector reads include true positives and false positives:
    detector(ii) = binornd(round(decays2),eff) + binornd(round(dt*noise),0.7);
    logdet(ii) = log(detector(ii));
end % End of For



% Save the results: Change the directory and file names as appropriate
%csvwrite("C:\Users\niked\Documents\Uni\2022S2\PVB303\Prac 2\Simulation\Times.csv", transpose(times));
%csvwrite("C:\Users\niked\Documents\Uni\2022S2\PVB303\Prac 2\Simulation\Detector.csv", transpose(detector));
%csvwrite("C:\Users\niked\Documents\Uni\2022S2\PVB303\Prac 2\Simulation\Detector_Log.csv", transpose(logdet));









