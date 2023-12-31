clear all
close all
clc

format long


% - - - - - - - - - - - - - - - - - - - 
% Model Parameters
% - - - - - - - - - - - - - - - - - - - 
Nx = 100;                       % number of x points
x = linspace(0,2*pi,Nx);        % x grid points
dx = x(2) - x(1);               % determine grid spacing
L = x(end) - x(1);              % determine the period
h = 0.5;                        % depth of the fluid
g = 1;                          % gravity

% - - - - - - - - - - - - - - - - - - - 
% Fourier Parameters
% - - - - - - - - - - - - - - - - - - - 
Nn = 32;                        % number of Fourier Modes
N = (-Nn:Nn)';                  % a vector of modes



% - - - - - - - - - - - - - - - - - - - 
% Initial Guess for Newton Method
% - - - - - - - - - - - - - - - - - - - 
iniGuessEta = .001*cos(x);
iniGuessEtaHat = hat(iniGuessEta,x,N);
iniGuessC = sqrt(g*tanh(h));
iniGuessX = [iniGuessEtaHat;iniGuessC];


% - - - - - - - - - - - - - - - - - - - 
% Use the internal Matlab "Newton's Method
% - - - - - - - - - - - - - - - - - - - 
options = optimset('TolFun',1e-8,'TolX',1e-8,'jacobian','on','display','iter');
X = fsolve(@(X) solveFunct(X,x,N), iniGuessX,options);


% - - - - - - - - - - - - - - - - - - - 
% Retreive and plot results
% - - - - - - - - - - - - - - - - - - - 
etaHat = X(1:end-1);
eta = invHat(etaHat,x,N);
plot(x,iniGuessEta,'.-')
hold on
plot(x,eta)
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% - - - - - - - - - - - - - - - - - - - 
% Compute the Fourier Coefficients
% - - - - - - - - - - - - - - - - - - - 
function uHat = hat(u,x,N)
dx = x(2)-x(1);
L = x(end)-x(1);
intTerm = exp(-1i*N*x).*u;      % create a matrix mesh of the Fourier integrand
uHat = sum(intTerm(:,2:end),2)*dx/L';   % This is the trapz. rule for periodic functions
end

% - - - - - - - - - - - - - - - - - - - 
% Compute the Fourier Series from Fourier Coefficients
% - - - - - - - - - - - - - - - - - - - 
function u = invHat(uHat,x,N)
sumTerm = exp(1i*N*x).*uHat;
u = sum(sumTerm,1);
if norm(imag(u))<1e-10
    u = real(u);
else
    disp('Complex function')
end
end


% - - - - - - - - - - - - - - - - - - - 
% Compute Derivative Spectrally
% - - - - - - - - - - - - - - - - - - - 
function uX = d(u,x,N)
uHat = hat(u,x,N);
uXHat = uHat.*(1i*N);
uX = invHat(uXHat,x,N);
end

% - - - - - - - - - - - - - - - - - - - 
% Compute the Nonlinear Eqns & Jacobian
% - - - - - - - - - - - - - - - - - - - 
function [F dF] = solveFunct(X,x,N);
Nn = max(N);
etaHat = X(1:end-1);
c = X(end);
g = 1;
h = 0.5;



eta = invHat(etaHat,x,N);
etaX = d(eta,x,N);

sqrtTerm = sqrt((c^2 - 2*g*eta).*(1 + etaX.^2));
F = zeros(2*Nn+2,1);
for n = -Nn:Nn
    
    if n==0
        F(Nn+1) = etaHat(Nn+1);
    else
        F(Nn+n+1) = hat(sqrtTerm.*sinh(n*(eta+h)),x,n);
    end
end

F(2*Nn+2) = sum(etaHat)-.001;
dF = jacobian(X,x,N);
end


% - - - - - - - - - - - - - - - - - - - 
% Calculate the Jacobian
% - - - - - - - - - - - - - - - - - - - 
function dF = jacobian(X,x,N)

Nn = max(N);
dF = zeros(2*Nn+2,2*Nn+2);

etaHat = X(1:end-1);
c = X(end);
g = 1;
h = 0.5;

eta = invHat(etaHat,x,N);
etaX = d(eta,x,N);
sqrtTerm = sqrt((c^2 - 2*g*eta).*(1 + etaX.^2));


for n = -Nn:Nn
    if n==0
        dF(Nn+1,Nn+1) = 1;
    else
        for k = -Nn:Nn
            sinhTerm = sinh(n*(eta+h));
            dFdEta = (sinhTerm./sqrtTerm.*(-g.*(1+etaX.^2)+(c^2 - 2*g*eta).*etaX.*1i*k)...
                + cosh(n*(eta+h)).*sqrtTerm.*n).*exp(1i*k*x);
            dF(n+Nn+1,k+Nn+1) = hat(dFdEta,x,n);
        end
    end
    dF(Nn+1,2*Nn+2) = hat(c./sqrtTerm.*(1+etaX.^2).*sinhTerm,x,n);
end
dF(2*Nn+2,1:2*Nn+1) = ones(1,2*Nn+1);

end