function res = dEke(obj)

% returns decomposed h1 energy matrix where split valence basis

% have to be executed after obj.solveHF

if(isempty(obj.Ehf)) % hartree fock not solved yet
    disp('please solveHF first');
    res = 0;
    return;
end
if(obj.singlezeta)
    res = obj.density .* obj.ke() .* obj.kemodmat;
else
    density_mat = obj.density();
    ke_mat = obj.ke();
    num_minbas = size(obj.minbas, 2);
    res = zeros(num_minbas,num_minbas);
    for i=1:num_minbas
        for j=1:num_minbas
            res(i,j) = sum(sum(density_mat(obj.minbas{i},obj.minbas{j}).*ke_mat(obj.minbas{i},obj.minbas{j})));
        end
    end
end
end
