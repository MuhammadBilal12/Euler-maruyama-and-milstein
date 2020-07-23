% Excercise 1
%Authors

% Muhammed Bilal 1736928



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% part (a)
% implementing euler and milstein on the SDE and plotting respective paths
rng('default')
mu = 0.06;          % Mu
sigma = 0.3;        % Volatility
S0 = 50;            % Initial Stock price
K = 90;             % Strike Price
N = 10000;          % Number of sample paths
M = 500;            % Number of discretization points
T = 1;              % Time period
dt = T/M;           % Step size
dW = randn(M,N);    % random numbers
delta = [0:dt:1];  % discretized time intervals
Seuler = zeros(M,N); % euler maruyama stock prices
Smilstein=zeros(M,N); % milstein stock rices
Seuler(1,:) = S0;       
Smilstein(1,:)=S0;
increment = 200;            % sample increment
%implementing the SDE using euler and milstein method over the discretized
%point
for i = 1:M         
    Seuler(i+1,:) = Seuler(i,:) + mu.*Seuler(i,:).*dt + sigma.*Seuler(i,:).*dW(i,:).*sqrt(dt); % Euler Maruyama
    %Milstein
    Smilstein(i+1,:) = Smilstein(i,:) + mu.*Smilstein(i,:).*dt - 0.5.*Smilstein(i,:)*((sigma)^2).*dt+ sigma.*Smilstein(i,:).*dW(i,:).*sqrt(dt) +0.5.*Smilstein(i,:).*((sigma)^2.*(dW(i,:).*sqrt(dt)).^2); %Milstein Scheme
end
%Plotting paths for euler maruyama by taking value of stock prices at some discretization point
 
figure;
hold on
plot(delta, Seuler(:,1), 'r-',delta, Seuler(:,2), 'r-',delta, Seuler(:,3), 'r-',delta,Seuler(:,4),'r-',delta,Seuler(:,5), 'r-', delta,Seuler(:,6),'r-')
legend('Sample Paths for Euler M')

xlabel('Time')
ylabel('Stock Price')
title('dS = uSdt + sigmaSdW implemented using Euler maruyama ')
hold off

%Plotting paths for milstein by taking value of stock prices at some discretization point
figure;
hold on
plot(delta, Smilstein(:,1), 'y-',delta, Smilstein(:,2), 'y-',delta, Smilstein(:,3), 'y-',delta,Smilstein(:,4),'y-',delta,Smilstein(:,5), 'y-', delta,Smilstein(:,6),'y-')
legend('Sample Paths for Milstein')
xlabel('Time')
ylabel('Stock Price')
title('dS = uSdt + sigmaSdW implemented using Milstein')
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%part(C)
% Calculate Option Prices using Euler and Milstein and comparing it with
% Black Scholes



Optioneuler = zeros(N/increment,1); %  option price for each sample set of euler
Optionmilstein = zeros(N/increment,1); %  option price for each sample set of milstein
Paths1 = zeros(N/increment,1); %  samples  of euler
Paths2 = zeros(N/increment,1); %   samples of milstein
for j = increment:increment:N
    eulerstock = Seuler(end,1:j)'; %  first i samples euler stock price
    milsteinstock=Smilstein(end,1:j)'; % first i samples milstein stock price
    Optioneuler(j/increment) = mean(max(eulerstock - K,0)); % Average Call Option price.
     Optionmilstein(j/increment) = mean(max(milsteinstock - K,0)); % Average Call Option price .
    Paths1(j/increment) = j; % Corresponding Number of samples of euler
    Paths2(j/increment) = j; % Corresponding Number of samples of milstein
end

% Black-Scholes Formula option price
d1 = (log(S0/K) + (mu + 0.5*sigma^2)*T)/(sigma*sqrt(T)); %conditional probability
d2 = d1 - sigma*sqrt(T); %probability that stock price is greater than strike price
N1 = 0.5*(1+erf(d1/sqrt(2))); % Commulative Normal Distribution
N2 = 0.5*(1+erf(d2/sqrt(2))); % Commulative Normal Distribution
C = S0*N1-K*exp(-mu*T)*N2; % Final price using Black Scholes
g=length(Paths1);
C1=zeros(g+1,1);
C1(:,:)=C;


%plotting the results of euler ,milstein and comparing with black scholes
figure;
plot(Paths1,Optioneuler,Paths2,Optionmilstein,[0; Paths1],C1)
hold on
plot(0,C, '-o')
ylim([C-0.1 C+0.2]);
xlabel('Number of Sample Paths N')
ylabel('Call Option Price')
legend('Euler Maruyama Price','Milstein Price','Black Scholes Price')
title('Comparison of Option Price between Euler-M,Milstein and BS')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% part(b)
% Computes the error of SDE between analytical solution,euler and milstein

rng('default')
clear dW;
N = 100; %  Sample Paths
M = 2000; %  Discretization Points
inc = 100;  % Increment in Discretization Points
SEuler = zeros(N,1); % stock price by Euler-M
SMilstein = zeros(N,1); % stock price by Milstein
STrue = zeros(N,1);  % stock price by Analytical Formula
Points = [inc:inc:M]'; %  number of Discretization Points
S0 = 50; %initial stock price
T = 1; %time period
r = 0.06; % interest rate
sigma = 0.3; %volatility
SE = S0; % Initial stock price for calculation for Euler-M
SM=S0;% Initial stock price for calculation for Milstein
ST = S0; % Initial stock price for calculation for Analytic
MeanErrorEuler = zeros(length(Points),1); %  Error between euler and analytic solution
MeanErrorMilstein = zeros(length(Points),1); % Error between milstein and analytic solution
dt = zeros(length(Points),1); % Vector which stores different step sizes used
for j = 1:length(Points)
    dt(j) = T/Points(j);  % Calculate step size for this particular number of discretization points
    SMilstein(:,:) = 0;
    SEuler(:,:) = 0;
    STrue(:,:) = 0;
    for i = 1:N
        SE = S0;
        SM=S0;
        ST = S0;
        dW = randn(Points(j),1)*sqrt(dt(j));
        W = cumsum(dW);
        for k = 1:Points(j)
            SE = SE + r*SE*dt(j) + sigma*SE*dW(k,1); % Euler Maruyama
             SM=SM + r*SM*dt(j)- 0.5*SM*((sigma)^2)*dt(j) +sigma*SM*dW(k,1) + 0.5*SM* ((sigma)^2)*(dW(k,1)^2); %Milstein
           
        end
            ST = S0*exp(((r-0.5*sigma^2)*T)+sigma*W(k)); % Analytic Solution
            SEuler(i) = SE;
            SMilstein(i)=SM;
           
            STrue(i) = ST;
    end
   
    MeanErrorEuler(j) = mean(abs(SEuler - STrue));
     MeanErrorMilstein(j) = mean(abs(SMilstein - STrue));
    SEuler(:,:) = 0;
    SMilstein(:,:)=0;
    STrue(:,:) = 0;
end

% plotting error of euler and milstein

figure;
ylim([0 0.3])
hold on
xlabel('Step Size')
ylabel('Mean Absolute Error')
title('Mean Abs Error of EulerMaruyama and Milstein from Analytical Solution ')

hold on
plot(dt,MeanErrorMilstein,'y-',dt,MeanErrorEuler,'r-')
legend('Euler Maruyama Error','Milstein Error')