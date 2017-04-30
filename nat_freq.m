
function [V_sorted,D_sorted]=nat_freq(M,K,R_TRA)

pai_G2TRA=[R_TRA,zeros(3);zeros(3),R_TRA];
M_prime=pai_G2TRA*M*pai_G2TRA';
K_prime=pai_G2TRA*K*pai_G2TRA';
A=inv(M_prime)*K_prime;
%Obtain eigenvalues and eigenvectors of A
[V,D]=eig(A);

[D_sorted, ind] = sort(diag(D),'ascend');
V_sorted = V(:,ind);

end