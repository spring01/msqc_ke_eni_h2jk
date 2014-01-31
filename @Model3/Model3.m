classdef Model3 < handle
    %{
   Interpolate between narrow and diffuse STO-3G matrix elements.
    %}
    properties
        % Input to the class
        frag;  % Fragment with regular STO-3G basis set
        
        % Most recent predictions of the model
        Ehf     % Hartree Fock energy
        Eorb    % (nbasis,1)      molecular orbital energies
        orb     % (nbasis,nbasis) molecular orbital coefficients

        Eelec   % electronic energy
        singlezeta
        minbas  % cell array projecting split valence basis onto minimal basis
        h2jk
        kemodmat
        enimod3mat
        h2jkmodnsqmat
        dEkestore
        dEenistore
        dEh2jkstore
        
        % Useful properties initialized in the constructor
        natom   % number of atoms in fragment
        nelec   % number of electrons in the fragment
        Z       % (1,natom) atomic numbers of the molecules
        aType   % (1,natom) atom type (initialized to Z)
        rcart   % (3,natom) cartesian coordinates of the atoms
        nbasis  % number of atomic (and molecular) basis functions
        
        %end
        %properties (Transient)
        densitySave   % most recent density matrix
        % used to start HF iterations
    end
    methods (Static)
        h2 = H2slater(F0, G1, F2)
    end
    methods
        
        function res = Model3(frag_)
            if (nargin ~= 0)
                res.frag = frag_;
                res.natom = frag_.natom;
                res.nelec = frag_.nelec;
                res.Z     = frag_.Z;
                res.aType = res.Z;
                res.rcart = frag_.rcart;
                res.nbasis = frag_.nbasis;
                
                res.singlezeta = res.initialize_singlezeta();
                res.genminbas();
                res.h2jk = reshape(res.h2(),res.nbasis.^2,res.nbasis.^2) - reshape(permute(res.h2(),[1 4 2 3]),res.nbasis.^2,res.nbasis.^2)./2;
                res.kemodmat = 1;
                res.enimod3mat = 1;
                res.h2jkmodnsqmat = 1;
            end
        end
        
        function res = initialize_singlezeta(obj)
            if( strcmpi(obj.frag.config.basisSet(1:3), 'sto') )
                res = 1;
            elseif ( strcmpi(obj.frag.config.basisSet, '6-31g') || strcmpi(obj.frag.config.basisSet, '3-21g') )
                res = 0;
            else
                res = 0;
            end
        end
        
        function genminbas(obj)
            % generate obj.minbas
            if( strcmpi(obj.frag.config.basisSet(1:3), 'sto') ) % minimal basis
                curr_ll_bas_start = 1;
                for atom_z=obj.frag.Z
                    if(atom_z>1) % heavy atom
                        obj.minbas{curr_ll_bas_start} = curr_ll_bas_start; % 1s
                        obj.minbas{curr_ll_bas_start+1} = curr_ll_bas_start+1; % 2s
                        obj.minbas{curr_ll_bas_start+2} = curr_ll_bas_start+2; % 2p1
                        obj.minbas{curr_ll_bas_start+3} = curr_ll_bas_start+3; % 2p2
                        obj.minbas{curr_ll_bas_start+4} = curr_ll_bas_start+4; % 2p3
                        curr_ll_bas_start = curr_ll_bas_start + 5; % ll jump to next atom
                    else % hydrogen
                        obj.minbas{curr_ll_bas_start} = curr_ll_bas_start; % 1s
                        curr_ll_bas_start = curr_ll_bas_start + 1; % ll jump to next atom
                    end
                end
            elseif( strcmpi(obj.frag.config.basisSet, '6-31g') || strcmpi(obj.frag.config.basisSet, '3-21g') ) % split valence
                curr_ll_bas_start = 1;
                curr_hl_bas_start = 1;
                for atom_z=obj.frag.Z
                    if(atom_z>1) % heavy atom
                        obj.minbas{curr_ll_bas_start} = curr_hl_bas_start; % 1s
                        obj.minbas{curr_ll_bas_start+1} = [curr_hl_bas_start+1 curr_hl_bas_start+5]; % 2s
                        obj.minbas{curr_ll_bas_start+2} = [curr_hl_bas_start+2 curr_hl_bas_start+6]; % 2p1
                        obj.minbas{curr_ll_bas_start+3} = [curr_hl_bas_start+3 curr_hl_bas_start+7]; % 2p2
                        obj.minbas{curr_ll_bas_start+4} = [curr_hl_bas_start+4 curr_hl_bas_start+8]; % 2p3
                        curr_ll_bas_start = curr_ll_bas_start + 5; % ll jump to next atom
                        curr_hl_bas_start = curr_hl_bas_start + 9; % hl jump to next atom
                    else % hydrogen
                        obj.minbas{curr_ll_bas_start} = [curr_hl_bas_start curr_hl_bas_start+1]; % 1s
                        curr_ll_bas_start = curr_ll_bas_start + 1; % ll jump to next atom
                        curr_hl_bas_start = curr_hl_bas_start + 2; % hl jump to next atom
                    end
                end
            end
        end
        
        function res = h1(obj)
            res = obj.ke() + obj.en();
        end
        
        function res = en(obj)
            res = sum(obj.frag.H1en,3);
        end
        
        function res = eni(obj)
            res = obj.frag.H1en;
        end
        
        function res = ke(obj)  % modifying matrix
            res   = obj.frag.KE;
        end
        
        function res = Hnuc(obj)
            res = obj.frag.Hnuc;
        end
        
        function res = h2(obj)
            res = obj.frag.H2;
        end
        
        function res = s(obj)
            res = obj.frag.S;
        end
        
    end % methods
    
end %
