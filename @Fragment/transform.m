function transform(obj,U)
% Unitary transformation, via matrix U(oldBasis,newBasis)

% Transform each property that changes with basis set
% H1(nbasis,nbasis) full H1 operator of fragment
obj.H1 = U' * (obj.H1* U);
% H1en(nbasis,nbasis,natom) electron-nuclear interaction
for iatom = 1:obj.natom
   obj.H1en(:,:,iatom) = U' * (squeeze(obj.H1en(:,:,iatom))*U);
end
% KE(nbasis,nbasis) kinetic energy
obj.KE = U' * (obj.KE * U);
% S(nbasis,nbasis) overlap
obj.S = U' * (obj.S * U);

% H2 (nbasis,nbasis,nbasis,nbasis) 2-elec interactions
N = obj.nbasis;
n2 = zeros(N,N,N,N);
for i=1:N
   for j=1:N
      for k=1:N
         for l=1:N
            for a=1:N
               for b=1:N
                  for c=1:N
                     for d=1:N
                        n2(i,j,k,l) = n2(i,j,k,l) + ...
                           obj.H2(a,b,c,d) * U(a,i) * U(b,j) ...
                           * U(c,k) * U(d,l);
                     end
                  end
               end
            end
         end
      end
   end
end
obj.H2 = n2;

% orb(nbasis,nbasis) molecular orbital coefficients
% orb is (atomic,molecular) so
obj.orb = U' * obj.orb;




