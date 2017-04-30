
function [V_sorted,D_sorted]=nat_freq(M,K,R_TRA)
%An example of Programming in MATLAB to obtain
%natural frequencies and mode shapes of MDOF
%systems
%Define [M] and [K] matrices
% M=[11 0;0 22]
% K=[1000 -500;-500 2000]
%Form the system matrix
pai_G2TRA=[R_TRA,zeros(3);zeros(3),R_TRA];
M_prime=pai_G2TRA*M*pai_G2TRA';
K_prime=pai_G2TRA*K*pai_G2TRA';
A=inv(M_prime)*K_prime;
%Obtain eigenvalues and eigenvectors of A
[V,D]=eig(A);
%V and D above are matrices.
%V-matrix gives the eigenvectors and
%the diagonal of D-matrix gives the eigenvalues
% Sort eigen-values and eigen-vectors

[D_sorted, ind] = sort(diag(D),'ascend');
V_sorted = V(:,ind);
%Obtain natural frequencies and mode shapes
% nat_freq_1 = sqrt(D_sorted(1))
% nat_freq_2 = sqrt(D_sorted(2))
% mode_shape_1 = V_sorted(:,1)
% mode_shape_2 = V_sorted(:,2)
end