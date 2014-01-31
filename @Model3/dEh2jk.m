function res = dEh2jk(obj)

% returns decomposed h2 energy quad-matrix where split valence basis
% have to be executed after obj.solveHF

if(isempty(obj.Ehf)) % hartree fock not solved yet
    disp('please solveHF first');
    res = 0;
    return;
end

if(obj.singlezeta)
    res = obj.density2pnsq() .* obj.h2jk() .* obj.h2jkmodnsqmat;
else
    % only for high level
    num_minbas = size(obj.minbas, 2);
    density2p_4mat = reshape(obj.density2pnsq(), obj.nbasis, obj.nbasis, obj.nbasis, obj.nbasis);
    h2jk_4mat = reshape(obj.h2jk, obj.nbasis, obj.nbasis, obj.nbasis, obj.nbasis);
    res = zeros(num_minbas,num_minbas,num_minbas,num_minbas);
    for i=1:num_minbas
        for j=1:num_minbas
            for k=1:num_minbas
                for l=1:num_minbas
                    res(i,j,k,l) = sum(sum(sum(sum(density2p_4mat(obj.minbas{i},obj.minbas{j},obj.minbas{k},obj.minbas{l}).*h2jk_4mat(obj.minbas{i},obj.minbas{j},obj.minbas{k},obj.minbas{l})))));
                end
            end
        end
    end
    res = reshape(res, num_minbas.^2, num_minbas.^2);
end
end

