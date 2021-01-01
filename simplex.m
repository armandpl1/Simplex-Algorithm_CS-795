function [x,optbasis,result]=simplex(A,b,c)
m=size(A,1)
n=size(A,2)
% Check correctness of input data(similar to simplex1 method
if size(b,1) ~= m || size(b,2) ~= 1
error('b wrong size')
end
if size(c,1) ~= n || size(c,2) ~= 1
error('c wrong size')
end
%checking b does not have negative entries and changing those negative to
%positive and flipping the sign of corresponding row in A
for i=1:m
if b(i,1)<0
b(i,1)=-1*b(i,1)
A(i,:)=-1*A(i,:)
end
end
%ADDITIONAL COEFFICIENTS FOR C TO CONSIDER XA
c_Xa=ones(m,1)
c_X=zeros(n,1)
c_phase1=[c_X;c_Xa]
%additional coefficient matrix to form phase 1 problem in SEF format
A_Xa=eye(m)
A_phase1=horzcat(A,A_Xa)
%initial BFS by assigfning Xa=b
init_BFS_phase1=[c_X;b]
initbasis_phase1=[n+1:m+n]
%phase 1 initialization
[phase1_x,phase1_basis,phase1_flag] = simplex1(A_phase1,b,c_phase1,init_BFS_phase1,initbasis_phase1)
% if at the end of phase 1, Xa had an entry not equal to zero the problem i
% infeasible and code terminates
for y=n+1:n+m
if phase1_x(y,1)~=0
x="infeasible"
optbasis="infeasible"
result=1
return
end
end
%selecting original variables for phase 2
phase2_x=phase1_x(1:n,:)
phase2_basis=zeros(1,n)
for k=1:n
phase2_basis(1,k)=ismember(k,phase1_basis)
end
% extending basis
basis=extendtobasis(A,phase2_basis)
% detecting columns that are associated with basic variables(logical value
% 1)
phase2base=zeros(m,1)
index=1
for l=1:n
if basis(1,l)==1
phase2base(index,1)=l
index=index+1
end
end
% running phase 2
[x,optbasis,result]=simplex1(A,b,c,phase2_x,phase2base)
%detecting unboundedness
if result==1
x="unbounded"
optbasis="unbounded"
result=2
return
end
end