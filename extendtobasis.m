function wmask = extendtobasis(A,bmask)
%  W = extendtobasis(A,B)
% This routine takes as input an m-by-n matrix A of rank m and bmask, a logical
% row vector of length n that indicates linearly independent columns of A.
% It returns wmask, also a logical row vector of length n with precisely m 1's
% that indicate a set of basis columns of A and that contain the w columns.

[m,n] = size(A);
A1 = A;

% Modified gram-schmidt loop

bmasksize = sum(bmask);
assert(bmasksize <= m)
bcols = find(bmask);
wmask = logical(zeros(1,n));

for j = 1 : m
    if j <= bmasksize
        curcol = bcols(j);
    else
        norms = sum(A1.^2, 1);
        notw = find(~wmask);
        [scrap,curcol0] = max(norms(notw));
        curcol = notw(curcol0);
    end
    wmask(curcol) = 1;
    A1(:,curcol) = A1(:,curcol) / norm(A1(:,curcol));
    A1(:,~wmask) = A1(:,~wmask) - A1(:,curcol) * (A1(:,curcol)'*A1(:,~wmask));
end
