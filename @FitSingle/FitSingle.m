classdef FitSingle < handle
    properties
        ll  % low level Model3
        hl  % high level Model3
        natom % number of atoms
        n   % number of minimal basis functions
        nsq
        nsqsq
        natom_nsq
        cn2
        c_cn2_2
        num_kepar
        num_enipar
        num_h2jkpar_a
        num_h2jkpar_b
        num_h2jkpar_c
        num_h2jkpar_d
        num_h2jkpar
        numpars
        parvec
    end
    methods
        
        function res = FitSingle(fileName, read_index)
            if (nargin < 2)
                read_index = 1;
            end
            load(fileName,'LL','HL');
            res.ll = Model3(LL{read_index,1});
            res.hl = Model3(HL{read_index,1});
            res.natom = res.ll.natom;
            res.n = res.ll.nbasis;
            res.nsq = res.n .^ 2;
            res.nsqsq = res.nsq .^ 2;
            res.natom_nsq = res.natom .* res.nsq;
            res.cn2 = res.n.*(res.n-1)./2;
            res.c_cn2_2 = res.cn2.*(res.cn2-1)./2;
            res.num_kepar = res.n + res.cn2;
            res.num_enipar = res.num_kepar .* res.natom;
            res.num_h2jkpar_a = res.num_kepar; 
            res.num_h2jkpar_b = res.cn2 + res.c_cn2_2;
            res.num_h2jkpar_c = res.num_h2jkpar_b;
            res.num_h2jkpar_d = res.n .* res.cn2;
            res.num_h2jkpar = res.num_h2jkpar_a + res.num_h2jkpar_b + res.num_h2jkpar_c+ res.num_h2jkpar_d;
            res.numpars = res.num_kepar + res.num_enipar + res.num_h2jkpar;
            res.hl.solvehf();
            res.hl.dEkestore = res.hl.dEke();
            res.hl.dEenistore = res.hl.dEeni();
            res.hl.dEh2jkstore = res.hl.dEh2jk();

            % MATLAB derping.
            warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
        end
        
        function solvehf(obj)
%             tic
            obj.ll.solvehf();
            obj.hl.solvehf();
%             toc
        end
        
        function res = err(obj, parvec)
            
            % store parameters
            obj.parvec = parvec;
            %parvec(2:length(parvec)) = ones(1,length(parvec)-1);

            % parse the input vector parvec into 3 sub vectors
            ke_parvec = parvec(1 : obj.num_kepar);
            eni_parvec = parvec(obj.num_kepar+1 : obj.num_kepar+obj.num_enipar);
            h2jk_parvec = parvec(obj.num_kepar+obj.num_enipar+1 : obj.num_kepar+obj.num_enipar+obj.num_h2jkpar);

            % generate (and store) modification matrices 
            obj.ll.kemodmat = obj.genkemodmat(ke_parvec);
            obj.ll.enimod3mat = obj.genenimod3mat(eni_parvec);
            obj.ll.h2jkmodnsqmat = obj.genh2jkmodnsqmat(h2jk_parvec);

            % solvehf for low level
            obj.ll.solvehf();
            
            % compute error matrices 
            Eke_errmat = 2 .* ( obj.ll.dEke() - obj.hl.dEkestore );
            Eeni_err3mat = 2 .* ( obj.ll.dEeni() - obj.hl.dEenistore );
            Eh2jk_errnsqmat = ( obj.ll.dEh2jk() - obj.hl.dEh2jkstore );

            % reshape error matrices into a vector, in order to execute lsqnonlin
            res(1 : obj.nsq) = reshape(Eke_errmat, 1, obj.nsq);
            res(obj.nsq+1 : obj.nsq+obj.natom_nsq) = reshape(Eeni_err3mat, 1, obj.natom_nsq);
            res(obj.nsq+obj.natom_nsq+1:obj.nsq+obj.natom_nsq+obj.nsqsq) = reshape(Eh2jk_errnsqmat, 1, obj.nsqsq);
%             res( n_2+n_2_times_natom+n_4+n_4+1 ) = (obj.llmodel.Ehf - obj.hlmodel.Ehf);  % total energy
%             disp(sum(res.^2));
        end
        
        function res = genhermitian(~, vec, dim)
            % input: vector and dimension of the result hermitian matrix 
            c_dim_2 = dim.*(dim-1)./2;
            resdiag = diag(vec(1 : dim));
            resoffdiag = zeros(dim) ;
            resoffdiag(tril(true(dim),-1)) = vec(dim+1 : dim+c_dim_2);
            res = resoffdiag + resoffdiag.' + resdiag;
        end
        
        function res = genkemodmat(obj, keparvec)
            % generate kemodmat from keparvec
            res = obj.genhermitian(keparvec, obj.n);
%             diag_keparvec = keparvec(1 : obj.n);
%             offdiag_keparvec = keparvec(obj.n+1 : obj.num_kepar);
%             resdiag = diag(diag_keparvec);
%             resoffdiag = zeros(obj.n, obj.n) ;
%             resoffdiag(tril(true(obj.n),-1)) = offdiag_keparvec ;
%             res = resoffdiag + resoffdiag.' + resdiag;
        end
        
        function res = genenimod3mat(obj, eniparvec)
            res = ones(obj.n,obj.n,obj.natom);
            for i=1:obj.natom
                res(:,:,i) = obj.genhermitian(eniparvec( (i-1).*obj.num_kepar+1:i.*obj.num_kepar ), obj.n);
            end
        end
        
        function res = genh2jkmodnsqmat(obj, h2jkparvec)
            % parse h2jkparvec into block a b c d 
            block_a = obj.genhermitian(h2jkparvec(1 : obj.num_h2jkpar_a), obj.n);
            block_b = obj.genhermitian(h2jkparvec(obj.num_h2jkpar_a+1 : obj.num_h2jkpar_a+obj.num_h2jkpar_b), obj.cn2);
            block_c = obj.genhermitian(h2jkparvec(obj.num_h2jkpar_a+obj.num_h2jkpar_b+1 : obj.num_h2jkpar_a+obj.num_h2jkpar_b+obj.num_h2jkpar_c), obj.cn2);
            block_d = reshape(h2jkparvec(obj.num_h2jkpar_a+obj.num_h2jkpar_b+obj.num_h2jkpar_c+1 : obj.num_h2jkpar), obj.n, obj.cn2);
            
            % generate result nsq by nsq matrix
            res = zeros(obj.nsq);
            scale1 = 1:obj.n;
            scale2 = obj.n+1:obj.n+obj.cn2;
            scale3 = obj.n+obj.cn2+1:obj.n+obj.cn2+obj.cn2;
            res(scale1, scale1) = block_a;
            res(scale2, scale2) = block_b;
            res(scale3, scale3) = block_b;
            res(scale2, scale3) = block_c;
            res(scale3, scale2) = block_c; % should be block_c' but block_c is hermitian
            res(scale1, scale2) = block_d;
            res(scale1, scale3) = block_d;
            res(scale2, scale1) = block_d';
            res(scale3, scale1) = block_d';
            
            % reorder the nsq matrix
            ordervec = zeros(1,obj.nsq);
            for i=1:obj.nsq
                p = ceil(i./obj.n);
                q = mod(i,obj.n);
                if(q==0)
                    q = obj.n;
                end
                if(p==q)
                    ordervec(i) = p;
                elseif(p<q)
                    ordervec(i) = obj.n + obj.n.*(obj.n-1)./2 + (2.*obj.n-p).*(p-1)./2 + (q-p);
                else
                    ordervec(i) = obj.n + (2.*obj.n-q).*(q-1)./2 + (p-q);
                end
            end
            res = res(ordervec,ordervec);
        end
        
    end
end