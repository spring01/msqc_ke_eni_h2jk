function res = dEeni(obj)

% returns decomposed h1 energy matrix where split valence basis 
% have to be executed after obj.solveHF

if(isempty(obj.Ehf)) % hartree fock not solved yet
    disp('please solveHF first');
    res = 0;
    return;
end

if(obj.singlezeta)
    for n=1:obj.natom
        density3mat(:,:,n) = obj.density;
    end
    res = density3mat .* obj.eni() .* obj.enimod3mat;
else
    density_mat = obj.density;
    eni_3mat = obj.eni();
    num_minbas = size(obj.minbas, 2);
    res = zeros(num_minbas,num_minbas);
    for i=1:num_minbas
        for j=1:num_minbas
            for n=1:obj.natom
                res(i,j,n) = sum(sum(density_mat(obj.minbas{i},obj.minbas{j}).*eni_3mat(obj.minbas{i},obj.minbas{j},n)));
            end
        end
    end
end
end
