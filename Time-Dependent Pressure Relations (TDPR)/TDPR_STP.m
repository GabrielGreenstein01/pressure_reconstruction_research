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
    rho =1;                         % fluid density

    % - - - - - - - - - - - - - - - - - - - 
    % Fourier Parameters
    % - - - - - - - - - - - - - - - - - - - 
    Nn = 10;                        % number of Fourier Modes
    N = (-Nn:Nn)';                  % a vector of modes

    % - - - - - - - - - - - - - - - - - - - 
    % Plotting Options
    % - - - - - - - - - - - - - - - - - - - 
    axes_options = {'interpreter','latex','fontsize',12};

    % - - - - - - - - - - - - - - - - - - - 
    % Initial Guess for Newton Method
    % - - - - - - - - - - - - - - - - - - - 
    iniGuessEta = 0.0001*cos(x);
    iniGuessEtaHat = hat(iniGuessEta,x,N);
    iniGuessC = sqrt(g*tanh(h));
    iniGuessX = [iniGuessEtaHat;iniGuessC];

    % - - - - - - - - - - - - - - - - - - - 
    % Use the internal Matlab "Newton's Method"
    % - - - - - - - - - - - - - - - - - - - 
    options = optimset('TolFun',1e-11,'TolX',1e-11,'jacobian','on','display','iter');
    [X,fval,exitflag,output,Jacobian] = fsolve(@(X) solveFunct(X,x,N), iniGuessX,options);

    % - - - - - - - - - - - - - - - - - - - 
    % Retreive and plot results
    % - - - - - - - - - - - - - - - - - - - 
    etaHat = X(1:end-1);
    eta = invHat(etaHat,x,N);
    c = X(end);

    % - - - - - - - - - - - - - - - - - - - 
    % Calculate the fk
    % - - - - - - - - - - - - - - - - - - - 
    fk = FC(eta, c, h, x, Nn, g);
    fk(Nn+1) = hat(-g*eta,x,0) %Resolves zero mode issue with test function: 1/2(z^2-x^2)

    % - - - - - - - - - - - - - - - - - - - 
    % Added better plotting options
    % - - - - - - - - - - - - - - - - - - - 
    
    subplot(2,1,1)
    stem(N, real(fk))
    xlabel('$N$',axes_options{:})
    title('Fourier Coefficients of $\hat\phi_x(x,-h)$',axes_options{:})

    subplot(2,1,2)
    fkInv = invHat(fk, x, N);
    P = -rho*(fkInv - g*h);
    
    plot(x, P-g*h,'-')
    hold on
    plot(x,rho*g*eta)
    plot(x,eta)
    leg = legend('$p- gh$','$\rho g \eta$','$\eta$','location','SouthEast');
    set(leg,axes_options{:});
    axis tight
    xlabel('$x$',axes_options{:})

    % - - - - - - - - - - - - - - - - - - - 
    % Store Data
    % - - - - - - - - - - - - - - - - - - - 
    save TDPR.mat

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % - - - - - - - - - - - - - - - - - - - 
    % Calculate the Fourier Coefficients of phi_x(x,-h)
    % - - - - - - - - - - - - - - - - - - - 
    function fk = FC(eta, c, h, x, Nn, g)
    fk = zeros(2*Nn+1, 1);
        for k = -Nn:Nn
            int = (g/k)*sinh(k*(eta + h)) + (c^2-2*g*eta).*cosh(k*(eta + h));
            fk(Nn+k+1) = hat(int, x, k);
        end
    end

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
    % Compute the Nonlinear Eqns & Jacobian
    % - - - - - - - - - - - - - - - - - - - 
    function [F dF] = solveFunct(X,x,N);
    Nn = max(N);
    etaHat = X(1:end-1);
    c = X(end);
    g = 1;
    h = 0.5;
    
    eta = invHat(etaHat,x,N);

    F = zeros(2*Nn+2,1);
    for n = -Nn:Nn

        if n==0
            F(Nn+1) = etaHat(Nn+1);
        else
            F(Nn+n+1) = hat((g/n)*cosh(n*(eta+h))+(c^2-2*g*eta).*sinh(n*(eta+h)),x,n);
        end
    end

    F(2*Nn+2) = sum(etaHat)-0.0001;
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

    for n = -Nn:Nn
        if n==0
            dF(Nn+1,Nn+1) = 1;
        else
            for k = -Nn:Nn
                dFdEta = exp(1i*k*x).*(n*(c^2-2*g*eta).*cosh(n*(eta+h))-g*sinh(n*(eta+h)));
                dF(n+Nn+1,k+Nn+1) = hat(dFdEta,x,n)*(abs(n-k) <= Nn);
            end
        end
        dF(n+Nn+1,2*Nn+2) = hat(2*c.*sinh(n*(eta+h)),x,n);
    end
    dF(2*Nn+2,1:2*Nn+1) = ones(1,2*Nn+1);
    dF(2*Nn+2, 2*Nn+2) = 0;
end
