function X = logm_agm_mp(A)

    n = size(A, 1); 
    epsilon = myeps(class(A));

    vareps = sqrt(epsilon) / norm(A,'fro');

    X = -mp('pi')/2 * inv(agm_Legendre_Taylor(vareps*A, epsilon));
    const = log(4/vareps);

    X(1:n+1:n^2) = X(1:n+1:n^2) + const;

end

% subfunctions

function X = agm_Legendre_Taylor(A, epsilon)
% This implements the Legendre variant for the computation of the arithmetic
% geometric mean of a matrix A such that ||A||^2 <= u(unit roundoff). It will
% provide innacurate results if the previous condition is not met because it
% approximates the inverse (I+A)^-1 with I-A

    mmax = 21;

    I=eye(length(A));
    k=1;
    [S, iterDB(k)]=sqrtm_dbp_mod(A,1,epsilon); % Denman and Beavers iteration for sqrtm
    P=2*S*(I-A);
    Q=(I+A)/4*(I+P); 
    curr_error(1)=norm(I-P,'fro');
    maxiter = 40;
    
%     while (k <= maxiter && (error(k) > 1 || ~test_cond(error(k),mmax,mp.eps())))
    while (k <= maxiter && curr_error(k) > 4 * myeps(class(A)))
        k=k+1;
        [S, iterDB(k)]=sqrtm_dbp_mod(P,1,epsilon);
    %     normP=norm(P,'fro');
        P=2*S/(I+P);   % MRHSLS
        Q=Q*(I+P)/2;
        curr_error(k)=norm(I-P,'fro');
    end
    
    X = Q;
    
    return
    
    %disp('number of iterations in each square root'), iterDB
    % switches to Taylor polynomials.
    
    degree=[2 3 5 7 9 13 17 21 ];
    Bound_error=find_bounds(degree, epsilon, curr_error(k)); 
    Delta0=I - P;
    poldegree=[]; % initializes a list of degrees of polynomials used
    m=degree(max(find(Bound_error <= curr_error(k)))); % the degree of the first one
    % Computation of coefficients c(j) in the Taylor series for Delta_{k+1}:
    rising=1/8; % this is the rising factorial of 2
    c(3)=1/8;   % the coefficient of Delta_k^2
    for j=3:m   % expresses each coefficient in terms of the previous one 
        rising=rising*(j-1.5)/j;
        c(j+1)=c(j)/2+rising;
    end

    while (k <= maxiter && curr_error(k) > epsilon)
            Delta1=polyvalm(fliplr(c(1:m+1)),Delta0);  % matrix polynomial evaluation.
            X=Q*(I-Delta1/2);
            k=k+1;
            curr_error(k)=norm(Delta1,'fro');
            Delta0=Delta1; Q=X;
            poldegree=[poldegree m];
            m=degree(max(find(Bound_error <= curr_error(k))));    
    end
    %disp('Errors and degrees of polynomials used'), error, poldegree
end

function [X,k] = sqrtm_dbp_mod(A, scale, epsilon)
% This is a modification of sqrtm_dbp from MFToolbox that avoids overflow 
% and underflow in the computation of the scaling factor
% k is the number of iterations

    n = length(A);
    if nargin < 2, scale = 1; end

    tol = mft_tolerance(A, epsilon);
    X = A;
    M = A;
    maxit = 25;

    for k = 1:maxit

       if scale == 1 
           [~, U]=lu(M); n=length(M);          % this avoids overflow and
           g=prod(abs(diag(U)).^(-1/(2*n)));   % underflow that may occur in
           % g = (abs(det(M)))^(-1/(2*n));   used in original sqrtm_dbp
           X = g*X; M = g^2*M;
       end

       Xold = X; invM = inv(M);

       X = X*(eye(n) + invM)/2;
       M = 0.5*(eye(n) + (M + invM)/2);

       Mres = norm(M - eye(n),'fro');

       reldiff = norm(X - Xold,'fro')/norm(X,'fro');
       if reldiff < 1e-2, scale = 0; end  % Switch to no scaling.

       if Mres <= tol, return; end

    end
    error('Not converged after %2.0f iterations', maxit)
end

function bounds = find_bounds(ms, u, delta_max)
    n = length(ms);
    p_local = mp.Digits();
    mp.Digits(2*p_local);
    delta_last = delta_max;
    for i = n:-1:1
        m = ms(i);
        delta_left = 0;
        delta_right = delta_last;
        found = false;
        while ~found
            if (abs(delta_right - delta_left) < u)
                bounds(i) = delta_left;
                delta_last = bounds(i);
                found = true;
                break;
            end
            delta_mid = (delta_right + delta_left) / 2;
            if test_cond(delta_mid, m, u)
                delta_left = delta_mid;
            else
                delta_right = delta_mid;
            end
        end
    end
    mp.Digits(p_local);
end

function satisfied = test_cond(deltak, m, u)
    % Check if (7.11) < u for a polynomial of degree m.
    d = mp.Digits();
    mp.Digits(4*d);
    ell2 = 1:m;
    p1_den = [1, cumprod(ell2)];
    p1_num = [1, 2.^ell2];
    p2 = [1, cumprod([0:m-1]-1/2)];
    ds = -sum((p1_num./p1_den).*p2)/(2^m);
%     ds*(deltak^(m+1)/(1-deltak))
    satisfied = ds*(deltak^(m+1)/(1-deltak)) < u;
    mp.Digits(d);
end

function e = myeps(curr_class)
%COMP_EPS    Machine epsilon.
%   COMP_EPS(CLASS) computes the machine epsilon of class CLASS.
    if(strcmp(curr_class, 'mp'))
        if (verLessThan('advanpix', '4.3.3.12177'))
            e = mp.eps();
        else
            e = eps(curr_class);
        end
    else
        e = eps(curr_class);
    end
end

function tol = mft_tolerance(A, epsilon)
%MFT_TOLERANCE   Convergence tolerance for matrix iterations.
%   TOL = MFT_TOLERANCE(A) returns a convergence tolerance to use in
%   the matrix iterations in the Matrix Function Toolbox applied to the
%   matrix A.  All functions in the toolbox call this function to set
%   the convergence tolerance.

    tol = sqrt(length(A))*epsilon/2;
end