%{
syms y(n) z
assume(n>=0 & in(n,"integer"))
f = y(n) - (2* y(n-1))- dirac(n)
fZT = ztrans(f,n,z)
syms pZT
fZT = subs(fZT,ztrans(y(n),n,z),pZT)
pZT = solve(fZT,pZT)
pSol = iztrans(pZT,z,n);
pSol = simplify(pSol)
pSol = subs(pSol,[y(0) y(1)],[1 2])

%{
nValues = 1:10;
pSolValues = subs(pSol,n,nValues);
pSolValues = double(pSolValues);
pSolValues = real(pSolValues);
stem(nValues,pSolValues)
title("COVID-19")
xlabel("Years (n)")
ylabel("Population p(n)")
grid on
%}
%}
%{
syms y(n) z

assume(n>=0 & in(n,"integer"))
Ro =1
num= [1]
den = [1 -Ro]
H = tf(num,den,1)
pzplot(H)
syms pZT
fZT = subs(fZT,H,pZT)
pZT = solve(fZT,pZT)
pSol = iztrans(pZT,z,n);
pSol = simplify(pSol)
pSol = subs(pSol,[y(0) y(1)],[1 2])
%%
r0 = [1,2,3,4]
b = 1;
for i=1:length(r0)
    a = [1,-r0(i)];
    for n=0:4
        infections = r0(i)^n
    end
    figure
    zplane(b,a)
    int2str(i);
    title("R0 = ",i)
end
%}
%{
% Transfer Function
R0 = 2.5; % Example R0 value
numerator = [1, 0];
denominator = [1, -R0];

H1 = tf(numerator, denominator,1)

% Pole-Zero Plot
pzmap(H1);
%{
r0 = [1,2,3]
b = 1;
for i=1:length(r0)
    a = [1,-r0(i)];
    for n=0:3
        infections = r0(i)^n
    end
    figure
    
    zplane(b,a)
    int2str(i);
    title("R0 = ",i)
end
%}

%2
syms a s
F = (s/((s-R0)));
f = ilaplace(F)
t= linspace(0,20);
f1 = (5*exp((5*t)/2))/2 + dirac(t);
figure
plot(f1)
% Initialize variables
n = 20; % Number of days
y = zeros(1, n);
y(1) = 1; % Assuming initial infection on day 1

% Iterate to find daily infections
for k = 2:n
    y(k) = 1 + R0 * y(k - 1);
end

% Display the results
%disp(y);

%3
threshold = 1e6; % 1 million
days_to_reach_threshold = ceil(log(threshold) / log(R0)) + 1;

disp(['It takes approximately ', num2str(days_to_reach_threshold), ' days to reach 1 million new daily infections.']);

%4
% one-point method
Y = 34285612;
dt = 365*1.5;
R01 = log(Y)/(548)

%linear-regression


%5
% Initialize variables
n = 20; % Number of days
y = zeros(1, n);
y(1) = 1; % Assuming initial infection on day 1

% Calculate new daily infections
for k = 2:n
    y(k) = 1 + R0 * y(k - 1);
end

% Plot new daily infections
figure;
subplot(2, 1, 1);
plot(1:n, y, 'bo-');
xlabel('Day');
ylabel('New Daily Infections');
title('New Daily Infections with R0 = 2.5');

% Integrate to calculate total infections
total_infections = cumsum(y);

% Plot total infections
subplot(2, 1, 2);
plot(1:n, total_infections, 'ro-');
xlabel('Day');
ylabel('Total Infections');
title('Total Infections with R0 = 2.5');
%}
%%
%{
% Initialize variables
n_max = 100; % Number of days
M = 12; % given limit
x = zeros(1, n_max+1);
y = zeros(1, n_max+1);
numerator = [1]; % coefficients
x(1) = 1; % initial infection
syms s 
ak = [-0.1,-0.1,-0.1,-0.15,-0.2,-0.25,-0.42,-0.34,-0.26,-0.25,-0.15,-0.1,1]; 
H1 = tf(numerator, ak)
F =  (-1)/(0.1*(s^12) + 0.1*(s^11) + 0.1*(s^10) + 0.15*(s^9) + 0.2*(s^8) + 0.25*(s^7) + 0.42*(s^6) + 0.34*(s^5) + 0.26*(s^4) + 0.25*(s^3)+ 0.15*(s^2)+ 0.1*(s) - 1);                                                                                                                   
f = ilaplace(F)
pzmap(H1);
% Calculate new daily infections
for n = 1:n_max
    for k = 1:M
        if n - k > 0
            y(n + 1) = y(n + 1) + ak(k) * y(n - k + 1);
        end
    end
    y(n + 1) = 1 - y(n + 1);
end

% Plot the daily infections
figure;
plot(0:n_max, y);
title('Daily Infections Over Time');
xlabel('Day');
ylabel('Daily Infections');
grid on;
figure;

% Integrate to calculate total infections
total_infections = cumsum(y);

% Plot total infections
plot(0:n_max, total_infections, 'ro-');
xlabel('Day');
ylabel('Total Infections');
title('Total Infections with R0 = 2.5');

%2
%first order was exponential while multiorder was a linear graph
% it never reaches 1 million
%}
%%
% Define the coefficients and parameters
M = 12;
ak = [0.1, 0.15, 0.25, 0.26, 0.34, 0.42, 0.25, 0.2, 0.15, 0.1, 0.1, 0.1];
n_max = 100;
rho_values = [0.25, 0.50, 0.75];
total_infections = zeros(length(rho_values), 1);

% Initialize arrays for the input (Kronecker delta) and output (daily infections)
x = zeros(1, n_max + 1);
y = zeros(1, n_max + 1);

% Loop over different ρ values
for rho_idx = 1:length(rho_values)
    rho = rho_values(rho_idx);
    
    % Set the Kronecker delta at day 0
    x(1) = 1;
    
    % Apply the IIR filter with scaled coefficients
    for n = 1:n_max
        for k = 1:M
            if n - k > 0
                y(n + 1) = y(n + 1) + (1 - rho) * ak(k) * y(n - k + 1);
            end
        end
        y(n + 1) = 1 - y(n + 1);
    end
    
    % Calculate the total number of infections for n = 100 days
    total_infections(rho_idx) = sum(y);
    disp('Total Infections for Different ρ Values:');
    disp(total_infections);
    % Plot the daily infections for the current ρ value
    subplot(1, length(rho_values), rho_idx);
    plot(0:n_max, y);
    title(['ρ = ' num2str(rho)]);
    xlabel('Day');
    ylabel('Daily Infections');
    grid on;
end

% Display the total number of infections for each ρ value


%%
%{
% Define parameters
R0 = 1.15;  % Reproduction number
K = 1e6;   % Population size
n_max = 100;  % Number of days

% Initialize arrays for the logistic and first-order models
x_logistic = zeros(1, n_max + 1);
x_first_order = zeros(1, n_max + 1);

% Initialize parameters for derivative calculation
D1_filter = [1, -1];  % First derivative filter coefficients
D2_filter = [1, -2, 1];  % Second derivative filter coefficients

% Apply the logistic model
for n = 0:n_max
    x_logistic(n + 1) = K / (1 + (K * (R0 - 1) - R0) * R0^(-(n + 1)) / (R0 - 1));
end

% Apply the first-order model
M = 12;  % Order of the first-order filter (adjust as needed)
ak = ones(1, M);  % Coefficients for the first-order model
x_first_order(1) = 1;  % Initial condition

for n = 1:n_max
    for k = 1:M
        if n - k > 0
            x_first_order(n + 1) = x_first_order(n + 1) + ak(k) * x_first_order(n - k + 1);
        end
    end
    x_first_order(n + 1) = 1 - x_first_order(n + 1);
end

% Plot the results
figure;
subplot(2, 1, 1);
plot(0:n_max, x_logistic, 'b', 'LineWidth', 2);
hold on;
plot(0:n_max, x_first_order, 'r--', 'LineWidth', 2);
title('Total Infections: Logistic vs. First-Order Model');
legend('Logistic Model', 'First-Order Model');
xlabel('Day');
ylabel('Total Infections');
grid on;

% Calculate first derivative and second derivative
first_derivative = conv(x_logistic, D1_filter, 'valid');
second_derivative = conv(x_logistic, D2_filter, 'valid');

% Find the inflection point
[~, max_derivative_index] = max(first_derivative);
zero_crossing_index = find(diff(sign(second_derivative)) == 2, 1);

disp(['Inflection point (First Derivative Maximum): Day ' num2str(max_derivative_index)]);
disp(['Inflection point (Zero-Crossing of Second Derivative): Day ' num2str(zero_crossing_index)]);
%}