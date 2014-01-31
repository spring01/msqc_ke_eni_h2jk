function solvehf(obj,eps,maxIter,minIter)

% eps: tolerance for the HF convergence
if (nargin < 2)
    eps = 1.0e-10;
end
if (nargin < 3)
    maxIter = 50000;
end
if (nargin < 4)
    minIter = 6;
end

[obj.orb,obj.Eorb,obj.Ehf,obj.Eelec] = obj.hf(eps,maxIter,minIter);
% inidensity = zeros(obj.nbasis);
% if(size(obj.densitySave,1) == 0)
%     inidensity = obj.frag.density;
% else
%     inidensity = obj.densitySave;
% end
% [obj.orb,obj.Eorb,obj.Ehf,obj.gstore,~,obj.densitySave] = parallelHFbeforeDIIS(obj.H1,obj.h2jk,obj.transx,obj.h1modmat,obj.gmodmat,obj.Hnuc,obj.frag.nelec,inidensity,eps,maxIter,minIter);

end

