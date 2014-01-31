s=dbstatus('-completenames');
save('myBreakpoints.mat', 's');
clear classes;
load('myBreakpoints.mat');
dbstop(s);

drive = 'c:';
load(fullfile(drive,'working\matlab\msqc_eemat\fitting\',...
   'h2geoms.mat'));
dataroot = fullfile(drive,'working\matlab\msqc_eemat\fitting\');
if (~exist(dataroot,'dir'))
   mkdir(dataroot);
end
% create linear array of zmats, with desired number of geometries
zmatsIn = zmats;
nopt = size(zmats,1);
ngeom = size(zmats,2);
% ngeom = 10;
zmats = cell(nopt*ngeom,1);
ic = 0;
for iopt = 1:nopt
   for igeom = 1:ngeom
      ic = ic+1;
      zmats{ic} = zmatsIn{iopt,igeom};
   end
end


LL = cell(length(zmats),1);
for izmat = 1:length(zmats)
   config1 = Config();
   config1.opt = 0; 
   config1.zmat = zmats{izmat};
   config1.basisSet = 'sto-3g';

   config1.silent = 1;
   frag1 = Fragment();
   success = frag1.runFragment(dataroot, config1);
   disp(['LL frag ',num2str(izmat)]);
   LL{izmat,1} = frag1;
%    LL{izmat,2} = frag;
%    LL{izmat,3} = frag;
end

HL = cell(length(zmats),1);
for izmat = 1:length(zmats)
   config2 = Config();
   config2.opt = 0; % do no geometry optimization
   config2.zmat = zmats{izmat};
   config2.method = 'hf';
   config2.basisSet = '6-31g';
   frag2 = Fragment();
   success = frag2.runFragment(dataroot, config2);
%    fragc = FragmentCached(frag2);
   disp(['HL frag ',num2str(izmat)]);
   HL{izmat,1} = frag2;
end

save(fullfile(dataroot,'h2.mat'),'LL','HL');
disp('Done');