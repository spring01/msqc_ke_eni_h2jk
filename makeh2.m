s=dbstatus('-completenames');
save('myBreakpoints.mat', 's');
clear classes;
load('myBreakpoints.mat');
dbstop(s);
clear;

drive = 'c:';
dataroot = fullfile(drive,'working\matlab\msqc_eemat\fitting\');
if (~exist(dataroot,'dir'))
    mkdir(dataroot);
end

% First we will create optimized versions of the z matrices

zmat1 = ZMatrix;
zmat1.file_name_base = 'fh';

zmat1.make_atom(9,0,0,0);
zmat1.make_atom(1,1,0,0);
zmat1.pars.bond_pars = [1.0];
config1 = Config();
config1.opt = 1; % do a geometry optimization
config1.zmat = zmat1;
config1.method = 'hf';
config1.basisSet = '6-31g';
config1.useCache = false;
frag1 = Fragment();
success1 = frag1.runFragment(dataroot, config1);

zmats = {};

% changing h-h
ngeom = 50;
for igeom = 1:ngeom
    zcurr = frag1.config.zmat.deepCopy;
    zcurr.pars.bond_pars(1) = 0.5 + 0.05.*(igeom-1);
    zmats{end+1} = zcurr.deepCopy;
end



save(fullfile(dataroot,'fhgeoms.mat'),'zmats');