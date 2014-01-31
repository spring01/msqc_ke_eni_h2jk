% wrapper of parallelHF to solve convergence issue
% When HF does not converge, projection operator eps_trial*|i><i| is added
% to generate a guess. While eps goes to zero, the Hartree-Forck equatrion
% is recovered.
% Directly solve. If failed, try the trick.
% Here, a self-adaptive scheme is employed.
% First, find a initial eps_trial that makes HF converge
% Second, try to decrease eps_trial towards zero
%   (a) If converged, success
%   (b) If not converged, backtrace and decrease the step
%   (c) If maximum number of calling parallelHF is reached, fail
function [orb,Eorb,Ehf,gstore,densityChanges,next_guessDensity] = parallelHFbeforeDIIS(H1,H2jk,X,h1modmat,gmodmat,Enuc,Nelec,guessDensity,eps,maxIter,minIter)

% suppress right division warning
warning('off','MATLAB:nearlySingularMatrix');
% default argument
if (nargin < 9)
    eps = 1.0e-8;
end
if (nargin < 10)
    maxIter = 1000;
end
if (nargin < 11)
    minIter = 5;
end

% direct solve
[orb,Eorb,Ehf,gstore,densityChanges,next_guessDensity] = ...
    parallelHF_e(H1,H2jk,X,h1modmat,gmodmat,Enuc,Nelec,guessDensity,eps,maxIter,minIter);

% if not converge
if (sum(densityChanges==0) == -1)
    current_eps_trial = 0.1; % this value does not matter the starting
    previous_eps_trial = 0.1; % previous eps, this is the starting
    current_step = 0.1; % starting step
    adaptive_factor = 2; % factor to decrease step
    maxtrial = 100; % maximum number of calling parallelHF
    iter_eps = 0;
    previous_Density = guessDensity;
    next_guessDensity = guessDensity;
    
    % search the first eps_trial
    isconverged = false;
    while(~isconverged && iter_eps <= maxtrial)
        iter_eps = iter_eps + 1;
        [orb,Eorb,Ehf,gstore,densityChanges,next_guessDensity] = ...
            parallelHF_e(H1,H2jk,X,h1modmat,gmodmat,Enuc,Nelec,previous_Density,eps,maxIter,minIter,previous_eps_trial,true);
        if (sum(densityChanges==0) == 0)
            % convergence failure
            previous_eps_trial = previous_eps_trial + current_step;
        else
            isconverged = true;
        end
    end
    
    disp(['parallelHF: init_eps_trial=', num2str(previous_eps_trial)]);
    
    previous_Density = next_guessDensity;
    isfinished = false;
    while(~isfinished && iter_eps <= maxtrial)
        iter_eps = iter_eps + 1;
        current_eps_trial = previous_eps_trial - current_step;
        if current_eps_trial < 0
            current_eps_trial = 0;
        end
        
        [orb,Eorb,Ehf,gstore,densityChanges,next_guessDensity] = ...
            parallelHF_e(H1,H2jk,X,h1modmat,gmodmat,Enuc,Nelec,previous_Density,eps,maxIter,minIter,current_eps_trial,true);
        if (sum(densityChanges==0) == 0)
            % not converge
            current_step = current_step / adaptive_factor;
            disp(['parallelHF: Guess not converge, current_eps_trial=', ...
                num2str(current_eps_trial), ' iter=', num2str(iter_eps), ' step=', num2str(current_step)]);
        else
            % converge
            previous_Density = next_guessDensity;
            previous_eps_trial = current_eps_trial;
            %
            if current_eps_trial == 0
                isfinished = true;
            end
            disp(['parallelHF: Guess converge, current_eps_trial=', ...
                num2str(current_eps_trial), ' iter=', num2str(iter_eps), ' step=', num2str(current_step)]);
        end
    end
    
    disp(['parallelHF: Terminated step=', num2str(current_step), ' iter=', num2str(iter_eps)]);
    
    if isfinished
        disp('parallelHF: Success with guess');
    end
    
    % detect maxtrial
    if iter_eps > maxtrial
        disp('parallelHF: FAILED, maxtrial is reached.');
        disp(['parallelHF: maxtrial=', num2str(maxtrial)]);
        disp(['parallelHF: current_eps_trial=', num2str(current_eps_trial)]);
    end
end
end

function [orb,Eorb,Ehf,gstore,densityChanges,P] = parallelHF_e(H1,H2jk,X,h1modmat,gmodmat,Enuc,Nelec,guessDensity,eps,maxIter,minIter,trial_eps)
% Solve Hartree Fock equations
% Input:
%   H1:    H1 matrix
%   H2:    H2 matrix
%   guessDensity:  starting density
%   eps:   convergence criteria
%          [defaults to 1e-8]
%   trial_eps: eps*|i><i|
% Output:
%   orb:   (nbasis,nbasis)
%            orbital coefficient matrix (atom basis #, mol orb #)
%   Eorb:  (nbasis,1)  molecular orbital energies
%   Ehf:    total Hartree-Fock energy

if (nargin < 9)
    eps = 1.0e-10;
end
if (nargin < 10)
    maxIter = 5000;
end
if (nargin < 11)
    minIter = 5;
end
if (nargin < 12)
    trial_eps = 0;
end

densityChanges = zeros(maxIter,1);

nbasis = size(H1,1);
% H2j {nbasis,nbasis} cell array of coulomb integrals
% H2k {nbasis,nbasis} cell array of exchange integrals
% H2j = cell(nbasis,nbasis);
% H2k = cell(nbasis,nbasis);
% for i=1:nbasis
%     for j=1:nbasis
%         H2j{i,j} = reshape(H2jsource(i,j,:,:), nbasis, nbasis);
%         H2k{i,j} = reshape(H2ksource(i,:,:,j), nbasis, nbasis);
%     end
% end


% Step 3 -- Calculate transformation matrix (eq. 3.167).
% X = transx;

% Step 4 -- Guess at density matrix -- all zeros right now.
Pn = guessDensity;
% Plast = Pn;  % previous density matrix
C = zeros(nbasis);
iter = 0;
finished = false;

%changeInDensityForCaution = sqrt(eps);
%changeInDensity = 10 * changeInDensityForCaution;

H1 = H1.*h1modmat;
% Begin iteration through.
while (~finished)  % Step 11 -- Test convergence
    
    %    if (changeInDensity > changeInDensityForCaution)
%     if issmix
%         % always use simple mixing
%         % simple mixing helps to converge in oscillating situation
%         P = 0.5 * Pn + 0.5 * Plast;
%     else
%         % simple mixing before maxIter/2
%         if (iter < maxIter/2)
%             P = 0.5 * Pn + 0.5 * Plast;
%         else
%             P = Pn;
%         end
%     end
    P = Pn;
    
    % update charges
    %    if (iter > 0)
    %       Q = zeros(size(obj.Z));
    %       GAP = zeros(1,obj.natom);
    %       P1 = P.*obj.S;
    %       GOP = sum(P1,1);
    %       for i = 1:obj.natom
    %          GAP(i) = sum(GOP(1,arange{i}));
    %          Q(i) = obj.Z(i)-GAP(i);
    %       end
    %       obj.charges(:,ienv+1) = Q;
    %       H1 = obj.H1(ienv);
    %    end
    
    % Step 5 -- Build 2-electron components of Fock matrix.
    G = reshape(H2jk*reshape(P,nbasis.*nbasis,1),nbasis,nbasis);
    
    % Step 5.5 Additional step for converging issue
    Prj_empty = zeros(nbasis);
    if trial_eps ~= 0
        unfilled = Nelec/2+1;
        % calc projection operator |i><i| for empty states
        for i = unfilled:nbasis
            Prj_empty = Prj_empty + C(:, i)*C(:, i)';
        end
        Prj_empty = trial_eps*Prj_empty;
    end
    
    % Step 6 -- Obtain F (fock matrix).
    F = H1 + G.*gmodmat + Prj_empty;
    
    % Step 7 -- Calculate the transformed F matrix.
    Ft = X'*F*X;
    
    % Step 8 -- Find e and the transformed expansion coefficient matrices.
    [Ct1,e1] = eig(Ft);
    e2 = diag(e1);
    [e, i1] = sort(e2);
    Ct = Ct1(:,i1);
    
    % Step 9 -- Transform Ct back to C.
    C = X*Ct; %#ok<MINV>
    
    % Step 10 -- Calculate the new density matrix.
    Plast = Pn;
    %Cj = conj(C);
    filled = 1:(Nelec/2);
    Pn = 2* C(:,filled)*( C(:,filled)');
    %     for i = 1:Nbasis
    %         for j = 1:Nbasis
    %             for a = 1:(Nelec/2)
    %                 Pn(i,j) = Pn(i,j) + (C(i,a)*Cj(j,a));
    %             end
    %             Pn(i,j) = Pn(i,j)*2;
    %         end
    %     end
    iter = iter + 1;
    %disp(['den change ',num2str( max(max(abs(P-Pn))))]);
    %disp(['diag Pn ',num2str(diag(Pn)')]);
    %Pn-Plast
    changeInDensity = max(max(abs(P - Pn)));
    densityChanges(iter) = norm(P-Pn);
    %disp(['iter ',num2str(iter),' del ',num2str(changeInDensity)]);
    if (iter > maxIter)
        finished = true;
    elseif (iter > minIter)
        if (changeInDensity < eps)
            finished = true;
        end
    end
end
% End of iteration of steps 5-11.

P = Pn;  % For convenience.

% Step 12: Output.

%Total energy
%3.184: E0 = 1/2 Sum(i,j) {P(j,i)[H1(i,j) + F(i,j)]}
%    Ee = 0;
%    for i = 1:Nbasis
%       for j = 1:Nbasis
%          Ee = Ee + P(j,i)*(H1(i,j)+F(i,j));
%       end
%    end
Ee = sum(sum(P.*(H1+F)));
Ehf = Ee/2 + Enuc;

% Orbital energies.
Eorb = e;

% Molecular orbital components.
orb = C;

gstore = G.*gmodmat;

% if (iter+1 > maxIter)
%     disp('You are living on the edge.. hartree fock didn''t converge');
% end
%{
Adapted from "Modern quantum chemistry", by Attila Szab?, Neil S. Ostlund
Numbered equations also adapted from here.
1. Specify a molecule
2. Calculate S(i,j), H^core (H1), and (i j|k l)(H2)
    -These first two steps are done by Gaussian
3. Diagonalize overlap matrix S and obtain X from 3.167
    3.167: X = S^(-1/2)
4. Guess the density matrix P (first guess is zeros here)
5. Calculate matrix G of 3.154 from P and H2
    G(i,j) = Sum(k, l){P(k,l)[(i j|l k)-1/2(i k|l j)]}
6. Add G to core-Hamiltonian  to get Fock matrix
    3.154: F(i,j) = H1(i,j) + G(i,j)
7. Calculate transformed Fock matrix F' = X'(t)FX
8. Diagonalize F' to obtain C' and epsilon
9. Calculate C = XC'
10. Form new density matrix P from C w/ 3.145
    3.145: P(i,j) = 2 Sum(1-Nelec/2){C(i,a) C*(j,a)}
11. Has P converged to within eps?
    No? -> Step 5 w/ new P from 10.
    Yes? -> Step 12
12. Use resultant solution, represented by C,P,F to calculate outputs
%}
end