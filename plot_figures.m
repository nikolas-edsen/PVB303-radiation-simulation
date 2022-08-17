%% Load data
load data1.mat;
load data2.mat;
load data3.mat;
load data4.mat;

%% Create detector noise plot
figure
set(gcf, "Position", [100,100,1100,400])
subplot(1,2,1)
plot(times1, detector1, '.')
ylim([0,15])
title('Counts Over Time Due to Detector Noise')
xlabel('Time (s)')
ylabel('Counts')

subplot(1,2,2)
histogram(detector1)
title('Distribution of Count Numbers Due to Detector Noise')
xlabel('Counts')
ylabel('Frequency')

text_x = get(gca, 'xlim') - 3;
text_y = get(gca, 'ylim') - 30;
text(text_x(2), text_y(2),"\sigma = " + std(detector1))



%% Create Secular Equilibrium Plot

figure
set(gcf, "Position", [200,300,1100,400])

subplot(1,2,1)
plot(times2, detector2, '.')
xlabel('Time (s)')
ylabel('Counts')
title('Counts Over Time at Secular Equilibrium')

subplot(1,2,2)
histogram(detector2)
title('Distribution of Count Numbers at Secular Equilibrium')
xlabel('Counts')
ylabel('Frequency')

text_x = get(gca, 'xlim') - 150;
text_y = get(gca, 'ylim') - 15;
text(text_x(2), text_y(2),"\sigma = " + std(detector2))



%% Create 137mBa decay plot

% Define curve fitting function
%F1 = @(p, t) p(1) * exp(-p(2) * t);

% Set initial parameters for least squares
%p0 = [7300, 0.004];

% Calculate best fit parameters for the given function
%[plsq, ~,~,~,~] = lsqcurvefit(F1, p0, times3, detector3);

% Plot graph of 137mBa decay data
figure

set(gcf, "Position", [200,300,1100,400])
subplot(1,2,1)
plot(times3, detector3, 'b.')
xlabel('Time (s)')
ylabel('Counts')
title('137mBa Decay')

% Find where the decay is below 2x the noise_level
noise_level = 9;

cutoff = 0;
for i = 1:length(detector3)
    if detector3(i) < (noise_level * 3)
        cutoff = i;
        break;
    end
end

t_lm = times3(1:cutoff);
y_lm = logdet3(1:cutoff);

% Create linear model
decay_model = fitlm(t_lm, y_lm)
c = decay_model.Coefficients{1,1};
c_e = decay_model.Coefficients{1,2};
m = decay_model.Coefficients{2,1};
m_e = decay_model.Coefficients{2,2};

decay_formula = @(t) c + m * t;
decay_str = sprintf('Formula: log(C) = %.4f t + %.4f', m, c)

% In log-linear coordinates with linear regression
subplot(1,2,2)
plot(times3, logdet3, 'b.')
hold on
plot(t_lm, decay_formula(t_lm), 'r-', 'LineWidth',1)
title('137mBa Decay (Logarithmic)')
xlabel('Time (s)')
ylabel('log(Counts)')
legend('Data', 'Linear Model')

text_x = get(gca, 'xlim') - 1230;
text_y = get(gca, 'ylim') - 1.5;
text(text_x(2),text_y(2),decay_str)



%% Create 137mBa Buildup plot
% Define curve fitting function
%F2 = @(p, t) p(1) * (1- exp(-p(2) * t));

% Set initial parameters for least squares
%p0 = [7300, 0.004];

% Calculate best fit parameters for the given function
%[plsq, ~,~,~,~] = lsqcurvefit(F2, p0, times4, detector4)

% Plot a graph of the 137mBa Buildup data
figure
set(gcf, "Position", [200,300,1100,400])

subplot(1,2,1)
plot(times4, detector4, 'b.')
xlabel('Time (s)')
ylabel('Counts')
title('137mBa Buildup')



tp = 30.17*525948.766*60; % T1/2 for 137Cs (in s)
lambda_p = log(2)/tp; % Decay rate constant for 137Cs (s^-1)

td = 2.552*60; % T1/2 for 137mBa (in s)
lambda_d = log(2)/td; % Decay rate constant for 137mBa (s^-1)

C_E = 5e14 * lambda_p/lambda_d * (1-exp(-lambda_d)) * 0.02; % Expected equilibrium count number from simulation
C_0 = 1e7 * 0.02 * (1-exp(-lambda_d)); % Expected initial count number

% Find cutoff for when count number exceeds the equilibrium
cutoff = 0;
for i = 1:length(detector4)
    if detector4(i) > C_E
        cutoff = i - 1;
        break;
    end
end


% logC = zeros(size(detector4));
% for i = 1:length(detector4)
%     if C_E > detector4(i)
%         logC(i) = log((C_E - detector4(i))/(C_E - C_0));
%     end
% end

logC = log((C_E - detector4)/(C_E - C_0));

t_lm = times4(1:cutoff);
y_lm = logC(1:cutoff);

% Create linear model
buildup_model = fitlm(t_lm, y_lm)
c = buildup_model.Coefficients{1,1};
c_e = buildup_model.Coefficients{1,2};
m = buildup_model.Coefficients{2,1};
m_e = buildup_model.Coefficients{2,2};

buildup_formula = @(t) c + m * t;
buildup_str = sprintf('Formula: !!log(C) = %.4f t + %.4f', m, c)



% In log-linear coordinates
subplot(1,2,2)
plot(times4, logC, 'b.');
hold on
plot(t_lm, buildup_formula(t_lm), 'r-', 'LineWidth',1)
title('137mBa Buildup (Logarithmic)')
xlabel('Time (s)')
ylabel('log(Counts)')
legend('Data', 'Linear Model')

text_x = get(gca, 'xlim') - 1230;
text_y = get(gca, 'ylim') - 2.5;
text(text_x(2),text_y(2),buildup_str)








