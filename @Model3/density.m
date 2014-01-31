function res = density(obj)
% Input:

% output:
%     density:  (nbasis,nbasis) 1-particle density matrix

%rho 1  from orb(atomic,mol)
%  psi_a = orb(i,a) phi_i
%  rho1(i,j) = 2 sum_(a occ) orb(i,a) orb(j,a)
%  rho1      = 2 orb(:,occ) * orb'(occ,j)

nocc = obj.nelec/2; %nocc: number of occupied orbitals

filledOrbs = obj.orb(:,1:nocc);

res = 2*filledOrbs*(filledOrbs'); % why is the 2 there??
% res = mm(2, filledOrbs);

end

