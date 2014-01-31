classdef Fragment < handle
    properties (SetAccess = private)
        config;        % structure holding calculation inputs
        dataPath;      % location of template and data files
        templateText;  % text from the template file
        gaussianFile;  % gaussian job (i.e. input) file (with charge keyword)
        fileprefix;    % prefix for files
        
        natom   % number of atoms in fragment
        nelec   % number of electrons in the fragment
        Z       % (1,natom) atomic numbers of the atoms
        rcart   % (3,natom) cartesian coordinates of the atoms
        npar    % number of parameters in template file
        
        nbasis  % number of atomic (and molecular) basis functions
        H1      % (nbasis,nbasis) full H1 operator of fragment
        H1en;   % (nbasis,nbasis,natom) electron-nuclear interaction
        KE;   % (nbasis,nbasis) kinetic energy
        H2;     % (nbasis,nbasis,nbasis,nbasis) 2-elec interactions
        S;      % (nbasis,nbasis) overlap
        Hnuc;   % nuclear-nuclear interaction energy
        
        Ehf     % Hartree Fock energy
        MP2     % MP2 Energy
        %CorrE   % Correlation Energy (MP2-Ehf)
        Eorb    % (nbasis,1)      molecular orbital energies
        orb     % (nbasis,nbasis) molecular orbital coefficients
        dipole  % (3,1)   dipole moment of molecule
        mulliken % (1,natom)  mulliken charge on the atoms

        basisAtom  % (nbasis,1) atom # on which the function is centered
        basisType  % (nbasis,1) l quantum number: 0=s 1=p 2=d 3=d etc
        basisSubType % (nbasis,1) m quantum number: s=1 p=3 d=6 (cartesian)
        basisNprims  % number of primitives in this function
        basisPrims   % {nbasis,1} cell array of matrices of size (2,nprims)
        %    with (1,:) being contraction coefficients and
        %         (2,:) being primimitive exponents
        
    end
    properties
        % TODO Need to do something about these
        gaussianPath = 'C:\ChemSoft\G09W';
        gaussianExe  = 'g09.exe';
    end
    methods (Access = private)
        initializeData(obj);
    end
    methods (Static)     
        [found,fileprefix] = findCalc(dataPath,config)
        [CorrE, MP2, Ehf, Eorb, orb, Nelectrons,  Z, rcart, ...
            dipole, mulliken, ...
            atom, type, subtype, nprims, prims ] = readfchk(fid1)
        [Eorb, orb, atom, Nelectrons, Ehf] = oldreadfchk(fid1)
        [S, H1, KE, H2, Enuc] = readpolyatom(fid1)
    end
    methods
        
        function dist = distance(obj, iatom, jatom)
            dist = norm(obj.rcart(:,iatom) - obj.rcart(:,jatom));
        end
        
    end % methods
end %
