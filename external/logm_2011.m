function [X, s_in, s, m] = logm_2011(A, maxsqrt)
%LOGM_NEW  Matrix logarithm by Schur-based inverse scaling and squaring.
%   X = LOGM_NEW(A,MAXSQRT) computes the logarithm of A, for a matrix
%   with no nonpositive real eigenvalues, using the inverse scaling and
%   squaring method with Pade approximation and a Schur decomposition.
%   [X, S, M] = LOGM_NEW(A) returns the number S of square roots
%   computed and the degree M of the Pade approximant.
%   At most MAXSQRT matrix square roots are  computed.
    
%   This code is intended for double precision.

%   Reference: A. H. Al-Mohy and N. J. Higham, Improved Inverse Scaling 
%   and Squaring Algorithms for the Matrix Logarithm, MIMS EPrint 2011.83,
%   The University of Manchester, October 2011.  
%   Name of corresponding algorithm in that paper: 
%   Algorithm 4.1/iss_schur_new.

%   Awad H. Al-Mohy and Nicholas J. Higham, October 19, 2011.

    if nargin < 2 || isempty(maxsqrt), maxsqrt = 100; end

    xvals = [1.586970738772063e-005
             2.313807884242979e-003
             1.938179313533253e-002
             6.209171588994762e-002
             1.276404810806775e-001
             2.060962623452836e-001
             2.879093714241194e-001];

    n = length(A);
    mmax = 7; foundm = false;
    % First form complex Schur form (if A not already upper triangular).
    if isequal(A,triu(A))
        T = A; Q = eye(n);
    else
        [Q,T] = schur(A,'complex');
    end
    if isdiag(T)      % Check if T is diagonal.
        d = diag(T);
        logd = log(d);
        X = (Q.*logd.')*Q';
        if isreal(logd)
            X = (X+X')/2;
        end
        s_in = 0;
        s = 0;
        m = 0;
    else
        T0 = T;
        if any( imag(diag(T)) == 0 & real(diag(T)) <= 0 )
            error('A must not have nonpositive real eigenvalues!')
        end
        p = 0;
        s0 = opt_cost(diag(T)); s = s0;
        for k = 1:min(s,maxsqrt)
            T = sqrtm_tri(T);
        end
        s_in = s0;
        
        d2 = normAm(T-eye(n),2)^(1/2); d3 = normAm(T-eye(n),3)^(1/3);
        alpha2 = max(d2,d3);
        if alpha2 <= xvals(2)
            m = find(alpha2 <= xvals(1:2),1); foundm = true;
        end
        
        while ~foundm
            more = 0; % more square roots
            if s > s0, d3 = normAm(T-eye(n),3)^(1/3); end
            d4 = normAm(T-eye(n),4)^(1/4);
            alpha3 = max(d3,d4);
            if alpha3 <= xvals(mmax)
                j = find( alpha3 <= xvals(3:mmax),1) + 2;
                if j <= 6
                    m = j; break
                else
                    if alpha3/2 <= xvals(5) && p < 2
                        more = 1; p = p+1;
                    end
                end
            end
            if ~more
                d5 = normAm(T-eye(n),5)^(1/5);
                alpha4 = max(d4,d5);
                eta = min(alpha3,alpha4);
                if eta <= xvals(mmax)
                    m = find(eta <= xvals(6:mmax),1) + 5;
                    break
                end
            end
            if s == maxsqrt, m = mmax; break, end
            T = sqrtm_tri(T); s = s + 1;
        end
        
        % Compute accurate superdiagonal of T^(1/2^s).
        for k = 1:n-1
            % Tkk = T0(k:k+1,k:k+1);
            % Tkk = powerm2by2(Tkk,1/2^s);
            % T(k:k+1,k:k+1) = Tkk;
            T(k:k+1,k:k+1) = powerm2by2(T0(k:k+1,k:k+1),1/2^s);
        end
        
        % Compute accurate diagonal of T^(1/2^s) - I.
        d = sqrt_power_1(diag(T0),s);
        % T = triu(T,1) + diag(d);
        T(1:n+1:end) = d;
        Y = logm_pf(T,m);
        
        X = 2^s*Y;
        
        % Compute accurate diagonal and superdiagonal of log(T).
        for l = 1:n-1
            X(l:l+1,l:l+1) = logm2by2(T0(l:l+1,l:l+1));
        end
        
        X = Q*X*Q';
    end
    if isreal(A), X = real(X); end
    
    % Nested functions

    %%%%%%%%%%%%%%%%%
    function s = opt_cost(d)
    % for i = 0:double(intmax)
    %     if max( abs(d-1)) <= xvals(mmax), s = i; return, end
    %     d = sqrt(d);
    % end
        s = 0;
        while norm(d-1,inf) > xvals(mmax)
            d = sqrt(d); s = s+1;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    function S = logm_pf(A,m)
    %LOGM_PF   Pade approximation to matrix log by partial fraction expansion.
    %   LOGM_PF(A,m) is an [m/m] Pade approximant to LOG(EYE(SIZE(A))+A).

        [nodes,wts] = gauss_legendre(m);
        % Convert from [-1,1] to [0,1].
        nodes = (nodes + 1)/2;
        wts = wts/2;

        S = zeros(n);
        for l1=1:m
            S = S + wts(l1)*(A/(eye(n) + nodes(l1)*A));
        end
    end

end

% Subfunctions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = sqrtm_tri(T)
% Compute upper triangular square root R of T, a column at a time.
    n = length(T);
    R = zeros(n);
    for j=1:n
        R(j,j) = sqrt(T(j,j));
        for i=j-1:-1:1
            R(i,j) = (T(i,j) - R(i,i+1:j-1)*R(i+1:j-1,j))/(R(i,i) + R(j,j));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = powerm2by2(A,p)
%POWERM2BY2    Power of 2-by-2 upper triangular matrix.
%   POWERM2BY2(A,p) is the pth power of the 2-by-2 upper
%   triangular matrix A, where p is an arbitrary real number.

    a1 = A(1,1);
    a2 = A(2,2);
    a1p = a1^p;
    a2p = a2^p;

    loga1 = log(a1);
    loga2 = log(a2);

    X = diag([a1p a2p]);

    if a1 == a2

        X(1,2) = p*A(1,2)*a1^(p-1);

    elseif abs(a1) < 0.5*abs(a2) || abs(a2) < 0.5*abs(a1)

        X(1,2) =  A(1,2) * (a2p - a1p) / (a2 - a1);

    else % Close eigenvalues.

        w = atanh((a2-a1)/(a2+a1)) + 1i*pi*unwinding(loga2-loga1);
        dd = 2 * exp(p*(loga1+loga2)/2) * sinh(p*w) / (a2-a1);
        X(1,2) = A(1,2)*dd;

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = sqrt_power_1(a,n)
%SQRT_POWER_1    Accurate computation of a^(2^n)-1.
%  SQRT_POWER_1(A,N) computes a^(2^n)-1 accurately.

% A. H. Al-Mohy.  A more accurate Briggs method for the logarithm.
% Numer. Algorithms, DOI: 10.1007/s11075-011-9496-z.

    if n == 0, r = a-1; return, end
    n0 = n;
    if angle(a) >= pi/2
        a = sqrt(a); n0 = n-1;
    end
    z0 = a - 1;
    a = sqrt(a);
    r = 1 + a;
    for i=1:n0-1
        a = sqrt(a);
        r = r.*(1+a);
    end
    r = z0./r;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = logm2by2(A)
%LOGM2BY2    Logarithm of 2-by-2 upper triangular matrix.
%   LOGM2BY2(A) is the logarithm of the 2-by-2 upper triangular matrix A.

    a1 = A(1,1);
    a2 = A(2,2);

    loga1 = log(a1);
    loga2 = log(a2);
    X = diag([loga1 loga2]);

    if a1 == a2
        X(1,2) = A(1,2)/a1;

    elseif abs(a1) < 0.5*abs(a2) || abs(a2) < 0.5*abs(a1)
        X(1,2) =  A(1,2) * (loga2 - loga1) / (a2 - a1);

    else % Close eigenvalues.
        dd = (2*atanh((a2-a1)/(a2+a1)) + 2*pi*1i*unwinding(loga2-loga1)) / (a2-a1);
        X(1,2) = A(1,2)*dd;

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,w] = gauss_legendre(n)
%GAUSS_LEGENDRE  Nodes and weights for Gauss-Legendre quadrature.
%   [X,W] = GAUSS_LEGENDRE(N) computes the nodes X and weights W
%   for N-point Gauss-Legendre quadrature.

% Reference:
% G. H. Golub and J. H. Welsch, Calculation of Gauss quadrature
% rules, Math. Comp., 23(106):221-230, 1969.

    i = 1:n-1;
    v = i./sqrt((2*i).^2-1);
    [V,D] = eig( diag(v,-1)+diag(v,1) );
    x = diag(D);
    w = 2*(V(1,:)'.^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = unwinding(z,k)
%UNWINDING    Unwinding number.
%   UNWINDING(Z,K) is the K'th derivative of the
%   unwinding number of the complex number Z.
%   Default: k = 0.

    if nargin == 1 || k == 0
        u = ceil( (imag(z) - pi)/(2*pi) );
    else
        u = 0;
    end
end