function [varargout] = logm_mp_abs(A, varargin)
%LOGM_MP  Multiprecision matrix logarithm and derivative
%   [X, S, M] = logm_mp(A,OPTS) computes the principal matrix logarithm of
%   the matrix A.
%
%   [Y, S, M] = logm_mp(A,E,OPTS) computes the Frechet derivative of the
%   matrix logarithm of A in the direction E.
%
%   [X, Y, S, M] = logm_mp(A,E,OPTS) computes both the principal matrix
%   logarithm of A and its Frechet derivative in the direction E.
%
%   [...] = logm_mp(...,'epsilon',EPSILON) specifies the tolerance to be
%   used to evaluate the Pade approximations. Default is 2^(-53).
%
%   [...] = logm_mp(...,'maxsqrtm', MAXSQRTM) specifies the maximum number
%   of square roots to be taken. Default is 100.
%
%   [...] = logm_mp(...,'maxdegree', MAXDEGREE) specifies the maximum
%   degree of the Pade approximants. Default is 100.
%
%   [...] = logm_mp(...,'approx', 'diagonal')
%   [...] = logm_mp(...,'approx', 'taylor') specifies the kind of Pade
%   approximant to use. Default is diagonal.
%
%   The output parameters S and M are the number of square roots taken by
%   the algorithm and the degree of the Pade approximation, respectively.

% The 'timing' option is an undocumented feature.

% Copyright (C) 2016 Massimiliano Fasi
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
% of the Software, and to permit persons to whom the Software is furnished to
% do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
% OTHER DEALINGS IN THE SOFTWARE.

    %% Parse and validate input.
    p = inputParser;
    addRequired(p, 'A', @ismatrix);
    addOptional(p, 'E', [], @(x)(ismatrix(x) && isnumeric(x)));
    addParameter(p,'epsilon', 2^-53, @(x)(x < 1));
    addParameter(p,'maxsqrt', 50, @(x)(x == round(x) && x > 0));
    addParameter(p,'maxdegree', 100, @(x)(x == round(x) && x > 0));
    addParameter(p,'approx', 'diagonal', @ischar);
    addParameter(p,'timing',false,@(x)(x == true || x == false));
    
    parse(p,A,varargin{:});

    A = p.Results.A;
    E = p.Results.E;
    epsilon = p.Results.epsilon;
    maxsqrt = p.Results.maxsqrt;
    maxdegree = p.Results.maxdegree;
    timing = p.Results.timing;

    if (strcmp(p.Results.approx, 'taylor'))
        usetaylor = true;
        eval_pade_bound = @scalar_bound_taylor;
        min_pade_deg = @min_taylor_deg;
        alpha_pade = @(A, m)(alpha(A,m,0));
    else
        usetaylor = false;
        eval_pade_bound = @scalar_bound_diagonal;
        min_pade_deg = @min_diagonal_deg;
        alpha_pade = @(A, m)(alpha(A,m,m));
    end
    [n1, n] = size(A);
    if (n ~= n1)
        error('The matrix ``A'' must be square.');
    end
    
    %% Set mp precision for single and double.
    if (isa(A, 'double'))
        mp.Digits(16);
    elseif (isa(A, 'single'))
        mp.Digits(8);
    end

    %% Parse and validate output.
    if (nargout > 4)
        error('This function returns at most four values.');
    end
    if (nargout < 4 && timing)
        error('When (''timing'', true) is used, the function returns four values.');
    end
    compute_frechet = false;
    if (~isempty(E))
        compute_frechet = true;
        [mE,nE] = size(E);
        if (nE ~= mE)
            error('The matrix ``E'' must be square.');
        end
        if (nE ~= n)
            error('The matrices ``A'' and ``E'' must have the same size.');
        end
    end
    compute_logm = ~(compute_frechet &&...
        (nargout == 3 || nargout == 1 || nargout == 0));

   %% Form complex Schur form (if A not already upper triangular).
    id_n = eye(n, class(A));
    schur_time = tic;
    if isequal(A, triu(A))
        T = A;
        Q = id_n;
    else
        d_old = mp.Digits();
%         mp.Digits(d_old + 0);
%         [Q, T] = schur(mp(A), 'complex');
%         mp.Digits(d_old);
%         Q = mp(Q);
%         T = mp(T);
        [Q, T] = schur(A, 'complex');
        if (compute_frechet)
            E = Q' * E * Q;
        end
    end
    T0 = T;
    time(1) = toc(schur_time);

    if any(imag(diag(T)) == 0 & real(diag(T)) <= 0)
        warning('``A'' must not have nonpositive real eigenvalues.')
    end

    %% Compute initial square roots.
    % Take square root until a pade approximant of degree maxdegree is
    % enough to achieve the required tolerance.
    init_sqrtm = tic;
    s = 0;
    Tm1 = T - id_n;
    a = -alpha_pade(Tm1, maxdegree);
    while ((max(abs(diag(Tm1))) >= 1 || abs(a) > 1 ||...
            eval_pade_bound(a, maxdegree) >= epsilon) &&...
            s < maxsqrt)
        [T, E, s] = compute_sqrtm(T, E, s, compute_frechet);
        Tm1 = T - id_n;
        a = -alpha_pade(Tm1, maxdegree);
    end
    time(2) = toc(init_sqrtm);

    %% Estimate good value for m.
    time(3) = 0;
    min_prob = tic;
    m = min_pade_deg(T, epsilon, maxdegree);
    
    old_m = m;
    count = 0;
    try_sqrtm = true;
    while (try_sqrtm && s < maxsqrt) % && count < 1)
        if usetaylor
            newm = max(ceil((sqrt(m) - 0.5)^2 - 1), 1);
        else
            newm = max(m - 2, 1);
        end
        a = -alpha_pade(Tm1, m) / 2;
        boundabs2 = eval_pade_bound(a, newm);

        % Check if taking sqrtm and a lower degree approximant will work.
        if boundabs2 > epsilon   % no
            try_sqrtm = false;
        else                                    % yes
            sub_sqrtm = tic;
            [T, E, s] = compute_sqrtm(T, E, s, compute_frechet);
            time(3) = time(3) + toc(sub_sqrtm);
            Tm1 = T - id_n;
            m = min_pade_deg(T, epsilon, m);
            if (m >= old_m)
                count = count + 1;
            else
                count = 0;
            end
            old_m = m;
        end
    end
    time(4) = toc(min_prob) - time(3);

    final_comp = tic;
    %% Compute accurate superdiagonal of T^(1/2^s).
    for l = 1:n-1
        T(l:l+1, l:l+1) = powerm2by2(T0(l:l+1, l:l+1), 1/2^s);
    end
    
    %% Compute accurate diagonal of T^(1/2^s) - I.
    T(1:n+1:end) = sqrt_power_1(diag(T0), s);
    
    if (~usetaylor || compute_frechet)
        [nodes, wts] = mp.GaussLegendre(m, 0, 1);
        nodes = cast(nodes, class(T));
        wts = cast(wts, class(T));
    end

    %% Compute matrix logarithm.
    if (compute_logm)
%         d_old = mp.Digits();
%         mp.Digits(d_old + 0);
%         T = mp(T);
        if (usetaylor)
            X = logm_taylor(T, m);
        else
            X = zeros(size(T,1), class(T));
            for k = 1:m
                X = X + wts(k) * ((id_n + nodes(k)*T) \ T);
            end
        end
        
        % DEBUG
%         old_d = mp.Digits();
%         mp.Digits(old_d * 2);
%         norm(X, 1)
%         norm(T)
%         norm(X - logm(mp(T) + eye(n)), 1) / norm(X, 1)
%         mp.Digits(old_d);

        X = 2^s * X;

        % Compute accurate diagonal and superdiagonal of log(T).
        for l = 1:n-1
            X(l:l+1, l:l+1) = logm2by2(T0(l:l+1, l:l+1));
        end
        
        X = Q * X * Q';
        if isreal(A)
            X = real(X);
        end
%         mp.Digits(d_old);
%         X = mp(X);
    end
    time(5) = toc(final_comp);

    %% Compute Frechet derivative.
    if (compute_frechet)
        L = zeros(n,n,class(T)); % The derivative
        for k = 1:m
            temp = id_n + nodes(k) * T;
            L = L + wts(k) * (temp \ E) / temp;
        end

        L = 2^s(1) * Q * L * Q';
        if (isreal(A))
            L = real(L);
        end
    end

    %% Prepare output.
    k = 1;
    if (compute_logm)
        [varargout{k}] = X;
        k = k + 1;
    end
    if (compute_frechet)
        [varargout{k}] = L;
        k = k + 1;
    end
    if (k < nargout)
        [varargout{k}] = s;
        k = k + 1;
        [varargout{k}] = m;
    end
    if (k < nargout)
        k = k + 1;
        [varargout{k}] = time;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [R, E, s] = compute_sqrtm(T, E, s, compute_frechet)
        if (ismp(T))
%             d_old = mp.Digits();
%             mp.Digits(d_old + 0);
            R = sqrtm_tri(mp(T));
%             mp.Digits(d_old);
        else
            R = sqrtm_tri_double(T);
        end
        if (compute_frechet)
            if (ismp(R))
                E = trisylv(R,R,E);
            else
                E = sylvester(R,R,E);
            end
        end
        s = s + 1;
    end

    function R = sqrtm_tri_double(T)
    %SQRTM_TRI Square root of quasi-upper triangular matrix.
    %   Works entirely in real arithmetic when possible.
        mm = size(T,1);
        switch mm
            case 1 % 1x1 block
                R = sqrt(T);
            case 2 % 2x2 block
                R = sqrtm_tbt(T);
            otherwise % Larger block - divide and conquer.
                n2 = floor(mm/2);
                if T(n2+1,n2) ~= 0
                    n2 = n2+1;
                end
                T11 = T(1:n2,1:n2);
                T22 = T(n2+1:end, n2+1:end);
                T12 = T(1:n2, n2+1:end);

                R11 = sqrtm_tri_double(T11);
                R22 = sqrtm_tri_double(T22);
                R12 = sylvester(R11, R22, T12);

                R(n2+1:mm, n2+1:mm) = R22;
                R(1:n2,1:n2) = R11;
                R(1:n2, n2+1:mm) = R12;
        end
    end

    function R = sqrtm_tbt(T)
    %SQRTM_TBT Square root of 2x2 matrix from block diagonal of Schur form.
    %
    %   Input T must be a 2x2 matrix from the block diagonal of a previously
    %   computed Schur form.

        if T(2,1) ~= 0
            % Compute square root of 2x2 quasitriangular block.
            % The algorithm assumes the special structure of real Schur form
            t11 = T(1,1); % t22 must equal to t11
            t12 = T(1,2);
            t21 = T(2,1);
            mu = sqrt(-t21*t12);
            if t11 > 0
                alpha = sqrt( (t11 + hypot(t11, mu))/2 );
            else
                alpha = mu / sqrt( 2*(-t11 + hypot(t11, mu)) );
            end
            R(2,2) = alpha;
            R(1,1) = alpha;
            R(2,1) = t21/(2*alpha);
            R(1,2) = t12/(2*alpha);
        else
            % Compute square root of 2x2 upper triangular block.
            t11 = T(1,1);
            r11 = sqrt(t11);
            t22 = T(2,2);
            r22 = sqrt(t22); 
            R(2,2) = r22;
            R(1,1) = r11;
            R(1,2) = T(1,2)/(r11 + r22);
        end
    end

    function X = powerm2by2(A, p)
    %POWERM2BY2    Power of 2-by-2 upper triangular matrix.
    %   POWERM2BY2(A, p) is the pth power of the 2-by-2 upper
    %   triangular matrix A, where p is an arbitrary real number.
        a1 = A(1, 1);
        a2 = A(2, 2);
        a1p = a1^p;
        a2p = a2^p;
        loga1 = log(a1);
        loga2 = log(a2);
        X = diag([a1p a2p]);
        if a1 == a2
            X(1,2) = p * A(1, 2) * a1^(p-1);
        elseif abs(a1) < 0.5 * abs(a2) || abs(a2) < 0.5 * abs(a1)
            X(1,2) =  A(1,2) * (a2p - a1p) / (a2 - a1);
        else % Close eigenvalues
            pi_typed = cast(mp('pi'), class(A));
            w = atanh((a2 - a1) / (a2 + a1)) +...
                1i * pi_typed * unwinding(loga2 - loga1, 0);
            X(1, 2) = 2 * A(1, 2) * exp(p * (loga1 + loga2) / 2) *...
                sinh(p * w) / (a2 - a1);
        end
    end

    function e = scalar_bound_diagonal(x, m)
    %LOG_DIAGONAL   Pade approximation to scalar logarithm.
    %   LOGM_PADE(A,M) computes the [M/M] Pade approximant to
    %   LOG(1+A) using bottom-up evaluation of the continued
    %   fraction expansion of the Pade approximants.
        d_old = mp.Digits();
%         mp.Digits(floor(d_old * 1.1));
        x = mp(x);
        y =  m / (4 * m - 2) * x;
        for j = m-1 : -1 : 1
            y =  (j * x) / ((j * 4 + 2) * (1 + y));
            y =  (j * x) / ((j * 4 - 2) * (1 + y));
        end
        y = x / (1 + y);
        e = abs(y - log(1 + x));
%         mp.Digits(d_old);
    end

    function S = logm_diagonal(A, m)
    %LOGM_DIAGONAL   Pade approximation to matrix logarithm.
    %   LOGM_PADE(A,M) computes the [M/M] Pade approximant to
    %   LOG(EYE(SIZE(A))+A) using the partial fraction expansion.
        [nodes, wts] = mp.GaussLegendre(m, 0, 1);
        nodes = cast(nodes, class(A));
        wts = cast(wts, class(A));
        id_n_loc = eye(size(A), class(A));
        S = zeros(size(A, 1), class(A));
        for j = 1:m
            S = S + wts(j) * (A / (id_n_loc + nodes(j) * A));
        end
    end

    function m = min_diagonal_deg(A, epsilon, maxd)
    %MIN_PADE_DEG   Degree of Pade approximation.
    %    MIN_PADE+DEG(T,EPSILON,MAXD) finds the minimum m such that the
    %    [m/m] Pade approximant is accurate enought to compute the
    %    matrix logarithm of (T + eye(n)) to accuracy EPSILON.
    
%         d_old = mp.Digits();
%         mp.Digits(floor(d_old));
        m_left = 1;
        m_right = maxd;
        found = false;
        AmI = A - eye(size(A,1), class(A));
        while (~found)
            m_middle = ceil((m_right + m_left) / 2);
            x = -alpha(AmI, m_middle, m_middle);

            if (abs(m_left - m_right) < 2)  
                m = m_right;
                found = true;
            elseif (x < -1 || scalar_bound_diagonal(x, m_middle) >= epsilon)
                m_left = m_middle;
            else % (abs(y - log (1 + x)) < epsilon)
                m_right = m_middle;
            end
        end
%         mp.Digits(d_old);
    end
    
    function e = scalar_bound_taylor(x, m)
    %LOGM_TAYLOR   Taylor approximation to matrix logarithm.
    %   LOGM_TAYLOR(A,M) computes the Taylor approximant of degree M to
    %   LOG(EYE(SIZE(A))+A) using Paterson-Stockmeyer algorithm.
        d_old = mp.Digits();
%         mp.Digits(floor(d_old * 1.1));
        x = mp(x);
        y = sum(-cumprod(repmat(-x, 1, m))./(1:m));
        e = abs(y - log(1 + x));
        mp.Digits(d_old);
    end

    function S = logm_taylor(A,m)
    %LOGM_TAYLOR   Taylor approximation to matrix logarithm.
    %   LOGM_TAYLOR(A,M) computes the Taylor approximant of degree M to
    %   LOG(EYE(SIZE(A))+A) using Paterson-Stockmeyer algorithm.

        v_0m = cast(0:1:m, class(A));
        c = (1./v_0m').*(-1).^(v_0m' + 1);
        c(1) = zeros(1, class(A));

        ss = ceil(sqrt(m));
        r = floor(m/ss);

        Xpow = cell(ss+1, 1);
        id_n_loc = eye(size(A,1), class(A));
        Xpow{1} = id_n_loc;
        Xpow{2} = A;
        for i=3:ss+1
            Xpow{i} = A*Xpow{i-1};
        end

        B = cell(r+1);
        for kk = 0 : r-1
            tt = c(ss*kk+1) * id_n_loc;
            for j=1:ss-1
                tt = tt + c(ss*kk+j+1)*Xpow{j+1};
            end
            B{kk+1} = tt;
        end

        B{r+1} = c(m+1)*Xpow{m-ss*r+1};
        for j=m-1:-1:ss*r
            if j == ss*r
               B{r+1} = B{r+1} + c(ss*r+1)*id_n_loc;
            else
               B{r+1} = B{r+1} + c(j+1)*Xpow{m-ss*r-(m-j)+1};
            end
        end

        As = Xpow{ss+1};
        S = zeros(size(A,1),class(A));
        for kk=r:-1:0
            S = S * As + B{kk+1};
        end
    end

    function m = min_taylor_deg(A, epsilon, maxd)
    %MIN_TAYLOR_DEG   Degree of Taylor approximation.
    %    MIN_TAYLOR_DEG(T,EPSILON,MAXD) finds the minimum M such that the
    %    Taylor series truncated to order M is accurate enought to compute
    %    the matrix logarithm of (T + eye(n)) to accuracy EPSILON.
        m_left = 1;
        m_right = maxd;
        found = false;
        AmI = A - eye(size(A), class(A));
        while (~found)
            m_middle = floor((m_right + m_left) / 2);
            x = -alpha(AmI, m_middle, 0);
            if (abs(m_left - m_right) < 2)  
                m = m_right;
                found = true;
            elseif (x < -1 || scalar_bound_taylor(x, m_middle) >= epsilon)
                m_left = m_middle;
            else % (abs(y - log (1 + x)) < epsilon)
                m_right = m_middle;
            end
        end
    end

    function X = logm2by2(A)
    %LOGM2BY2    Logarithm of 2-by-2 upper triangular matrix.
    %   LOGM2BY2(A) is the logarithm of the 2-by-2 upper triangular matrix
    %   A.
        a1 = A(1, 1);
        a2 = A(2, 2);
        loga1 = log(a1);
        loga2 = log(a2);
        X = diag([loga1 loga2]);
        if a1 == a2
            X(1, 2) = A(1, 2)/a1;
        elseif abs(a1) < 0.5 * abs(a2) || abs(a2) < 0.5 * abs(a1)
            X(1, 2) =  A(1, 2) * (loga2 - loga1) / (a2 - a1);
        else % Close eigenvalues.
            pi_typed = cast(mp('pi'), class(A));
            X(1, 2) = A(1, 2) * (2 * atanh((a2 - a1) / (a2 + a1)) +...
                2 * pi_typed * 1i * unwinding(loga2 - loga1, 0)) / (a2 - a1);
        end
    end

    function r = sqrt_power_1(a,n)
    %SQRT_POWER_1    Accurate computation of a^(2^n)-1.
    %  SQRT_POWER_1(A,N) computes a^(2^n)-1 accurately.

    % A. H. Al-Mohy.  A more accurate Briggs method for the logarithm.
    % Numer. Algorithms, DOI: 10.1007/s11075-011-9496-z.
        if (n == 0)
            r = a - 1;
        else
            n0 = n;
            if imag(a) < 0
                a = sqrt(a);
                n0 = n - 1;
            end
            z0 = a - 1;
            a = sqrt(a);
            r = 1 + a;
            for i = 1 : n0 - 1
                a = sqrt(a);
                r = r .* (1 + a);
            end
            r = z0 ./ r;
        end
    end

    function y = alpha(A, k, m)
    % ALPHA    Compute max(||X^p||^(1/p), ||X^(p+1)||^(1/(p+1)))
        R = double(A);
        d = floor((1 + sqrt(2 * (m + k) + 5)) / 2);
        dp1 = normAm(R, d)^(1/d);
        d = d + 1;
        dp2 = normAm(R, d)^(1/d);
        y = mp(max(dp1, dp2));
    end

    function u = unwinding(z,k)
    %UNWINDING    Unwinding number.
    %   UNWINDING(Z,K) is the K'th derivative of the
    %   unwinding number of the complex number Z.
    %   Default: k = 0.
        if k == 0
            pi_typed = cast(mp('pi'), class(A));
            u = ceil((imag(z) - pi_typed) / (2 * pi_typed));
        else
            u = 0;
        end
    end
end