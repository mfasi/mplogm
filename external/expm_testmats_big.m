function [A, nmats] = expm_testmats_big(k,n)
%EXPM_TESTMATS  "Big" test matrices for matrix exponential.
%   [A, NMATS] = EXPM_TESTMATS_BIG(K,N) selects the K'th test matrix.
%    NMATS is the number of test matrices available.
%    N sets (approximately) the dimension of those matrices for which the
%    dimension is variable.
%    These matrices are mainly discretizations of operators.
%    EIGTOOL is required.

if nargout > 1, nmats = 15; A = 0; end
if nargin < 1, return; end
if nargin < 2, n = 8; end

switch k

   case 1
   % EIGTOOL.  4n^2-by-4n^2
   A = convdiff_fd_demo(round(sqrt(n)/2));  % Nearest to n-by-n.

   case 2
   % EIGTOOL.
   A = transient_demo(n);

   case 3
   % EIGTOOL.
   A = davies_demo(n);

   case 4
   % EIGTOOL.
   A = basor_demo(n);

   case 5
   % EIGTOOL.
   A = airy_demo(n);

   case 6
   % EIGTOOL.
   A = chebspec_demo(n);

   case 7
   % EIGTOOL.
   A = companion_demo(n);

   case 8
   % EIGTOOL.
   A = convdiff_demo(n+1);

   case 9
   % EIGTOOL.
   A = demmel_demo(n);

   case 10
   % EIGTOOL.
   A = hatano_demo(n);

   case 11
   % EIGTOOL.
   A = landau_demo(n);

   case 12
   % EIGTOOL.  Dimension N^2/2 + 3*N/2 + 1;
   N = n; n = 1+round((-3 + sqrt(1+8*N))/2);
   A = markov_demo(n);

   case 13
   % EIGTOOL.
   [A,P] = riffle_demo(n);

   case 14
   % EIGTOOL.
   A = supg_demo(round(sqrt(n)));

   case 15
   % EIGTOOL.
   A = twisted_demo(n);

end
