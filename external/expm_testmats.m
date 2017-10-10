function [A, nmats] = expm_testmats(k,n)
%EXPM_TESTMATS  Test matrices for matrix exponential.
%   [A, NMATS] = EXPM_TESTMATS(K,N) selects the K'th test matrix.
%    NMATS is the number of test matrices available.
%    N sets the dimension of those matrices for which the dimension
%    is variable.

% Authors: Awad H. Al-Mohy and Nicholas J. Higham

if nargout > 1, nmats = 32; A = 0; end
if nargin < 1, return; end
if nargin < 2, n = 8; end

switch k

   case 1
   % \cite[Test 1]{ward77}.
   A = [4 2 0; 1 4 1; 1 1 4];

   case 2
   % \cite[Test 2]{ward77}.
   A = [29.87942128909879     .7815750847907159 -2.289519314033932
          .7815750847907159 25.72656945571064    8.680737820540137
        -2.289519314033932   8.680737820540137  34.39400925519054];

   case 3
   % \cite[Test 3]{ward77}.
   A = [-131 19 18;
        -390 56 54;
        -387 57 52];

   case 4
   % \cite[Test 4]{ward77}.
   A = gallery('forsythe',10,1e-10,0);

   case 5
   % \cite[p. 370]{naha95}.
   T = [1 10 100; 1 9 100; 1 11 99];
   A = T*[-0.001 0 0; 0 -1 0; 0 0 -100]/T;

   case 6
   % \cite[Ex.~2]{kela98}.
   A = [0.1 1e6; 0 0.1];

   case 7
   % \cite[p.~655]{kela98}.
   A = [0  3.8e3 0    0   0
        0 -3.8e3 1    0   0
        0 0     -1  5.5e6 0
        0 0      0 -5.5e6 2.7e7
        0 0      0   0   -2.7e7];

   case 8
   % \cite[Ex.~3.10]{dipa00}
   w = 1.3; x = 1e6; n = 8;
   A = (1/n) * [w*ones(n/2) x*ones(n/2)
                zeros(n/2)  -w*ones(n/2)];

   case 9
   A = rosser;
   A = 2.05*A/norm(A,1);  % Bad case for expm re. cost.

   case 10
   A = [0 1e4;
        -1e4 0];  % exp = [cos(x) sin(x); - sin(x) cos(x)], x = 100;

   case 11
   A = 1e2*triu(randn(n),1);  % Nilpotent.

   case 12 % log of Cholesky factor of Pascal matrix. See \cite{edst03}.
   A = zeros(n); A(n+1:n+1:n^2) = 1:n-1;

   case 13 % \cite[p.~206]{kela89}
   A = [48 -49 50 49; 0 -2 100 0; 0 -1 -2 1; -50 50 50 -52];

   case 14 % \cite[p.~7, Ex I]{pang85}
   A = [0    30 1   1  1  1
       -100   0 1   1  1  1
        0     0 0  -6  1  1
        0     0 500 0  1  1
        0     0 0   0  0  200
        0     0 0   0 -15 0];

   case 15 % \cite[p.~9, Ex II]{pang85}
   % My interpretation of their matrix for arbitrary n.
   % N = 31 corresponds to the matrix in above ref.
   A = gallery('triw',n,1);  m = (n-1)/2;
   A = A - diag(diag(A)) + diag(-m:m)*sqrt(-1);
   for i = 1:n-1, A(i,i+1) = -2*(n-1)-2 + 4*i; end

   case 16 % \cite[p.~10, Ex III]{pang85}
   A = gallery('triw',n,1,1);
   A = A - diag(diag(A)) + diag(-(n-1)/2:(n-1)/2);

   case 17
   % \cite[Ex.~5]{kela89}.
   A = [0 1e6; 0 0];   % Same as case 6 but with ei'val 0.1 -> 0.

   case 18
   % \cite[(52)]{jemc05}.
   g = [0.6 0.6 4.0]; b = [2.0 0.75];
   A = [-g(1)       0    g(1)*b(1)
          0        -g(2) g(2)*b(2)
        -g(1)*g(3)  g(3) -g(3)*(1-g(1)*b(1))];

   case 19
   % \cite[(55)]{jemc05}.
   g = [1.5 0.5 3.0 2.0 0.4 0.03]; b = [0.6 7.0];
   A1 = [-g(5)     0      0
          0      -g(1)    0
        g(4)     g(4)   -g(3)];
   A2 = [-g(6)       0    g(6)*b(2)
          0        -g(2)  g(2)*b(1)
          0         g(4) -g(4)];
   A = [zeros(3) eye(3); A2 A1];

   case 20
   % \cite[Ex.~3]{kela98}.
   A = [-1 1e7; 0 -1e7];

   case 21
   % \cite[(21)]{mopa03}.
   Thalf = [3.8235*60*24 3.10 26.8 19.9]/60;  % Half lives in seconds/
   a = log(2)./Thalf;  % decay constant
   A = diag(-a) + diag(a(1:end-1),-1);

   case 22
   % \cite[(26)]{mopa03}.
   a1 = 0.01145;
   a2 = 0.2270;
   A = [-a1              0  0
         0.3594*a1     -a2  0
         0.6406*a1     a2  0];

   case 23
   % \cite[Table 1]{kase99}.
   a = [4.916e-18
        3.329e-7
        8.983e-14
        2.852e-13
        1.373e-11
        2.098e-6
        9.850e-10
        1.601e-6
        5.796e-8
        0.000];
   A = diag(-a) + diag(a(1:end-1),-1);

   case 24
       % Jitse Niesen sent me this example.
       lambda = 1e6 * 1i;
       mu = 1/2*(-1+sqrt(1+4*lambda));
       A = [ 0, 1; lambda, -1 ] - mu*eye(2);

   case 25 % Awad

       A = [1 1e17;0 1];

   case 26 % Awad

       b = 1e3; x = 1e10;
       A = [ 1-b/2   b/2 ; -b/2   1+b/2 ];
       A = [A          x*ones(2);
            zeros(2)       -A    ];

   case 27 % Awad
       b = 1e5;
       A = [ 1-b/2   b/2 ; -b/2   1+b/2 ];

   case 28 % Awad
       b = 1e3;
       A = [ 1-b/2   b/2 ; -b^4/2   1+b/2 ];

   case 29
   % EIGTOOL.
   A = godunov_demo/10;

   case 30
   % \cite[(14.17), p. 141]{trem05}.
   A = 10*[0 1 2; -0.01 0 3; 0 0 0];

   case 31
   A = triu(schur(gallery('invol',13),'complex'),1);

   case 32
   % \cite{kuda10}
   alpha = 1; beta = 1;  % No values are given in the paper, unfortunately.
   A = -eye(n) + alpha/2*(diag(ones(n-1,1),1) + diag(ones(n-1,1),-1));
   A(1,2) = beta; A(n,n-1) = beta;

end
