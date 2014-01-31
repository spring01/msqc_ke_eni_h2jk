function res = density2pnsq(obj)

% For Restricted Hartree Fock theory, the 2-particle density is derivable 
% from the 1-particle density, as follows:

columned_dens = reshape(obj.density(), obj.nbasis.^2, 1);
res = columned_dens*columned_dens';


% res = reshape(res, obj.nbasis, obj.nbasis, obj.nbasis, obj.nbasis);
% 
% p1hf = obj.density;
% nb = obj.nbasis;
% p2hf = zeros(nb, nb, nb, nb);
% for a=1:nb
%    for b=1:nb
%       for c=1:nb
%          for d=1:nb
% %             p2hf(a,b,c,d) = 2.*p1hf(a,b).*p1hf(c,d) - p1hf(a,d).*p1hf(c,b);
%             p2hf(a,b,c,d) = p1hf(a,b).*p1hf(c,d);
%          end
%       end
%    end
% end
% 
% res = p2hf;

end


