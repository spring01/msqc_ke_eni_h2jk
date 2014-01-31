function bool = compareMolecule(obj,ztarget)
%If they are the same molecule the output will be 1 otherwise it will be 0
atoms = obj.atoms;
tatoms = ztarget.atoms;
natoms = length(atoms);
tnatoms = length(tatoms);
try
    bool = 0;
    if natoms ~=tnatoms;
        return
    end
    for iatom = 1:natoms
        if (atoms{iatom}.compare_atoms(tatoms{iatom}))
            bool = 1;
        else
            bool = 0;
            return
        end
    end
catch exception
    disp('error caught');
    exception.stack
end
end
