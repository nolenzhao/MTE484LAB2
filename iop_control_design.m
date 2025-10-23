
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This template is for the IOP control design with SPA, and is incomplete.
% You need to complete it by replacing every * with the correct code.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set time step
T = 0.01;

%% Plant Poles and Coefficients in its Partial Fraction Decomposition
% T = 2 * pi / ( 1/tau * ~7) 


K1 = 2.065 ;
tau = 0.0231;

% G = tf([K1], [tau, 1, 0]);
s = tf('s');
G = K1 / (s * (s * tau + 1));
Gd = c2d(G, T);
Gd_poles = pole(Gd);
Gd_zeros = zero(Gd);

%% Get Poles

stableRealPlantPoles = [0.6486];
stableComplexPlantPoles = [];
unstablePlantPoles = [1];
% 
stablePlantPoles = [stableRealPlantPoles stableComplexPlantPoles];
qs = [stablePlantPoles unstablePlantPoles];

% coefficents go in order of the poles
cs = [0.0166629 0.0205519];


n = length(qs);
nhat = length(stablePlantPoles);
nreal = length(stableRealPlantPoles);
ncomplex = length(stableComplexPlantPoles);

% verify that your plant is correct!
% z = tf('z', T)
% G = 0;
% for k=1:n
%     G = G + cs(k) / (z - qs(k));
% end
% G

%% Poles Chosen in the Simple Pole Approximation of W[z]

realWPoles = [];
complexWPoles = [-0.1645 - 0.0702i  -0.1645 + 0.0702i   0.0750 - 0.2416i   0.0750 + 0.2416i   0.3067 - 0.0442i   0.3067 + 0.0442i   0.2476 + 0.2583i   0.2476 - 0.2583i  -0.0290 + 0.3989i ...
                 -0.0290 - 0.3989i  -0.3231 + 0.2960i  -0.3231 - 0.2960i  -0.4728 + 0.0217i  -0.4728 - 0.0217i  -0.4171 - 0.2864i  -0.4171 + 0.2864i  -0.1895 - 0.5021i  -0.1895 + 0.5021i ...
                 0.1204 - 0.5527i   0.1204 + 0.5527i   0.4094 - 0.4294i   0.4094 + 0.4294i   0.5944 - 0.1751i   0.5944 + 0.1751i   0.6301 + 0.1376i   0.6301 - 0.1376i   0.5130 + 0.4299i ... 
                 0.5130 - 0.4299i   0.2748 + 0.6360i   0.2748 - 0.6360i  -0.0304 + 0.7149i  -0.0304 - 0.7149i  -0.3398 + 0.6546i  -0.3398 - 0.6546i  -0.5957 + 0.4702i  -0.5957 - 0.4702i ...
                 -0.7543 + 0.1975i  -0.7543 - 0.1975i  -0.7916 - 0.1159i  -0.7916 + 0.1159i];

ps = [realWPoles complexWPoles];

mreal = length(realWPoles);
mcomplex = length(complexWPoles);
m = length(ps);

%% Calculation of alpha, beta, gamma, and gamma hat

alpha = zeros(m);

for i=1:m
    for k=1:n
        alpha(i,i) = alpha(i,i) + cs(k) / (ps(i) - qs(k));
    end
end

beta = zeros(n,m);

for i=1:m
    for k=1:n
        beta(k,i) = cs(k) / (qs(k) - ps(i));
    end
end

gamma = zeros(n-nhat,m);

for i=1:m
    for j=(nhat+1):n
        gamma(j-nhat,i) =cs(j) / (qs(j) - ps(i));
    end
end

gammaHat = zeros(n-nhat,nhat);

for k=1:nhat
    for j=(nhat+1):n
        gammaHat(j-nhat,k) = cs(j) / (qs(j) - qs(k));
    end
end

% verify on a simple example that alpha, beta, gamma, and gammahat are correct!
alpha
beta
gamma
gammaHat

%% Determination of A and b matrices for IOP equations

A = [alpha eye(m) zeros(m,nhat);
     beta [zeros(nhat,m) eye(nhat);
           zeros(size(beta,1)-nhat,m+nhat)];
     zeros(size(gamma)) gamma gammaHat];

b = [zeros(m+size(beta,1),1);
     -cs(nhat+1:n)'];

%% Determination of step response matrices

% time horizon
K = 100;

step_ry = zeros(K,m+nhat);

for k=1:K
    for i=1:m
        step_ry(k,i) = -(1-ps(i)^k) / (1-ps(i));
    end
    for j=1:nhat
        step_ry(k,m+j) = -(1-qs(j)^k) / (1-qs(j));
    end
end

step_ru = zeros(K,m);

for k=1:K
    for i=1:m
        step_ru(k,i) = (1-ps(i)^k) / (1-ps(i));
    end
end

% verify on a simple example that step_ry and step_ru are correct!
step_ry
step_ru

%% Determination of steady state vector

steadyState = zeros(1,m+nhat);

for i=1:m
    steadyState(i) = 1/ (1-ps(i));
end

for k=1:nhat
    steadyState(m+k) = 1/(1-qs(k));
end

% verify on a simple example that steadyState is correct!
steadyState

%% Defining the variables for the optimization

wreal = sdpvar(mreal,1,'full');
wcomplex = sdpvar(mcomplex/2,1,'full','complex');
w = wreal;
for i=1:(mcomplex/2)
    w = [w;
         wcomplex(i);
         conj(wcomplex(i))];
end

xreal = sdpvar(mreal,1,'full');
xcomplex = sdpvar(mcomplex/2,1,'full','complex');
x = xreal;
for i=1:(mcomplex/2)
    x = [x;
         xcomplex(i);
         conj(xcomplex(i))];
end

xhatreal = sdpvar(nreal,1,'full');
xhatcomplex = sdpvar(ncomplex/2,1,'full','complex');
xhat = xhatreal;
for i=1:(ncomplex/2)
    xhat = [xhat;
            xhatcomplex(i);
            conj(xhatcomplex(i))];
end


%% Defining the objective function and constraints for the optimization

%Objective = 0;
Objective = 0;

% IOP constraint
Constraints = [A * [w;x;xhat] == b];

% input saturation constraint
Constraints = [Constraints,
               max(step_ru * w) <= 6, 
               min(step_ru * w) >= -6];

% steady state constraint
Constraints = [Constraints,
               1 + steadyState * [x;xhat] == 0];

% overshoot constraint
Constraints = [Constraints,
               max(step_ry * [x;xhat]) <= 1.05 * (-steadyState * [x;xhat])];

% settling time constraint
Ts = 0.25;
jhat = ceil(Ts/T);
Constraints = [Constraints,
               max(step_ry(jhat:end, :) * [x;xhat]) <= 1.02 * (-steadyState * [x;xhat]),
               min(step_ry(jhat:end, :) * [x;xhat]) >= 0.98 * (-steadyState * [x;xhat])];

%% Solving the optimization problem

% set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','mosek');

% solve the problem
sol = optimize(Constraints,Objective,options);

% obtain the solution
wsol = value(w);
xsol = value(x);
xhatsol = value(xhat);

%% Plotting the solution

figure(1)
plot(T*(1:K),step_ry*[xsol;xhatsol]);
xlabel('Time [s]');
ylabel('y[k]');

figure(2)
plot(T*(1:K),step_ru*wsol);
xlabel('Time [s]');
ylabel('u[k]');

% use log scale for heat map?
log_scale_flag = 1;

% heat map
figure(3)
t = linspace(0,2*pi);
plot(cos(t),sin(t),'k--');
hold on;
if log_scale_flag
    scatter(real(ps),imag(ps),50,log(abs(wsol)),'filled');
    scatter(real(qs(1:nhat)),imag(qs(1:nhat)),50,log(abs(xhatsol)),'filled');
else
    scatter(real(ps),imag(ps),50,abs(wsol),'filled');
    scatter(real(qs(1:nhat)),imag(qs(1:nhat)),50,abs(xhatsol),'filled');
end
hold off;
colormap(jet);
colorbar
%% Recover the transfer functions

z = tf('z',T);

% calculate W
W = 0;
for i=1:m
    W = W + wsol(i) / (z- ps(i));
end

% calculate X
X = 1;
for i=1:m
    X = X + xsol(i)/ (z - ps(i));
end
for k=1:nhat
    X = X + xhatsol(k) / (z - qs(k));
end

% remove the imaginary coefficients in W
[num,den] = tfdata(W);
num{1} = real(num{1});
den{1} = real(den{1});
W = tf(num,den,T);

% remove the imaginary coefficients in X
[num,den] = tfdata(X);
num{1} = real(num{1});
den{1} = real(den{1});
X = tf(num,den,T);

% find the poles and zeros of W and X
zpk(W)
zero(W)
pole(W)
zpk(X)
zero(X)
pole(X)

%% Verify design in discrete time

% compute D by hand
j = sqrt(-1);
D = ;

% compute T_ry and T_ru by hand
T_ry = *;
T_ru = *;

figure(1)
hold on;
*;
hold off;

figure(2)
hold on;
*;
hold off;

