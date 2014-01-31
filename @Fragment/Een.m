function res = Een(obj,iatom)
% energy of interaction of the molecule 


res = sum(sum( obj.density.*obj.H1en(:,:,iatom) ));
% res = elementWiseCombine(obj.density, obj.H1en(:, :, iatom));
end

