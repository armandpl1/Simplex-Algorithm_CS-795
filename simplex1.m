function [finalx, finalbfs, unboundedflag] = ...
  simplex1(A,b,c,initx,initbasis)
% [finalx, finalbfs, unboundedflag] = SEFlinprog1(A,b,c,initx,initbasis)
% This function solves a linear programming problem in standard equality
% form using the simplex method given an initial BFS.
%
% The input arguments are as follows.
% A: an m-by-n matrix of rank m (constraint matrix)
% b: an m-vector (RHS of constraints)
% c: an n-vector (cost vector).
% (initx, initbasis): an initial BFS; x is the n-vector of the BFS
% that must be feasible (i.e., x>=0 and A*x=b), and
% initbasis is a vector containing m integers lying in 1:n.  The
% matrix A(:,initbasis) must be invertible.
% 
% The return arguments are as follows.
% finalx: the optimal solution to the LP.
% finalbfs: the basis that goes with final x, a subset of m integers lying in
% 1:n.
% unboundedflag: set to 0 if the LP is bounded (i.e., an optimal found), or
% 1 if it is unbounded.

m = size(A,1);
n = size(A,2);

% Check correctness of input data
if size(b,1) ~= m || size(b,2) ~= 1
    error('b wrong size')
end
if size(c,1) ~= n || size(c,2) ~= 1
    error('c wrong size')
end

if size(initx,1) ~= n || size(initx,2) ~= 1
    error('initx wrong size')
end

if length(initbasis) ~= m || any(initbasis < 1) || any(initbasis > n)
    error('initbasis invalid')
end

if any(initx < 0)
    error('initial x infeasible for sign constraints')
end

if norm(A * initx - b) >= 1e-15 * norm(b)
    error('initial x infeasible for equality constraints')
end

basisflag = logical(zeros(n,1));
basisflag(initbasis) = 1;
nonbasiccols = find(basisflag == 0);
if length(nonbasiccols) ~= n - m
    error('initbasis has duplicate entries');
end
if any(initx(nonbasiccols) > 0)
    error('initx has nonzero entries in nonbasic columns')
end

x = initx;
B = initbasis;
N = nonbasiccols;
AB = A(:,B);
invAB = inv(AB);

while 1 % infinite loop until return statement encountered
    % Note parentheses in next statement: (c(B)' * invAB) * A requires
    % O(m*n) ops, but c(B)' * (invAB * A) requires O(m^2 * n) ops!
    reducedcosts = c' - (c(B)' * invAB) * A;  
    % The following test catches the case when an optimizer is found.
    if all(reducedcosts(N) >= 0)
        finalx = x;
        finalbfs = B;
        unboundedflag = 0;
        return
    end
    % Let pos be the position with respect to N of the minimum
    % reduced cost (Dantzig rule)
    [scrap, pos] = min(reducedcosts(N));
    % j is the column index with respect to A of the column entering the
    % basis
    j = N(pos);
    d = zeros(n,1);
    d(j) = 1;
    d(B) = -invAB * A(:,j);
    Bsub = find(d(B) < 0);
    if length(Bsub) == 0
        finalx = x;
        finalbfs = B;
        unboundedflag = 1;
        return
    end
    % B0 is the set of column indices from B with d(B) negative.
    B0 = B(Bsub);
    [theta,pos2] = min(-x(B0) ./ d(B0));
    % pos2 is the position of the min with respect to B0.
    i = B0(pos2);
    % Now i is the variable leaving the basis.
    x = x + theta * d;
    % Let pos3 be the position of i with respect to B.
    pos3 = find(B == i);
    % Change this entry of B from i to j; update N accordingly.
    B(pos3) = j;
    N(pos) = i;
    % Update AB.
    AB(:,pos3) = A(:,j);
    % Update invAB using the Sherman-Morrison-Woodbury formula.
    u = A(:,j) - A(:,i);
    v = zeros(m,1);
    v(pos3) = 1;
    % Note parentheses in SMW formula to ensure efficient computation.
    invAB = invAB - (invAB * u) * (v' * invAB) / (1 + v' * invAB * u);
end
end
