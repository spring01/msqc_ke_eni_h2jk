function bool = initializeZipData( fragment,zipFileName, envTarget)
    % runs gaussian in scratch directories, and puts all relavant files
    % into zipFileName
    if (nargin < 3)
       envTarget = [];
    end
    
%  SET VARS  
    
    jobname = 'full';
    gjf_file = [jobname,'.gjf'];
    
    gaussianPath = fragment.gaussianPath;

%  FILE MANAGEMENT / RUN GAUS
    if (~fragment.config.silent)
       disp('initializing the data');
    end    
    % Do the calculation
    if (isempty(envTarget))
       gjf_text = buildGjf( fragment );
    else
       gjf_text = buildGjf( fragment, 0, fragment.config.charge, ...
          fragment.config.spin, {envTarget.gaussianKeywords} );
       extraText = envTarget.gaussianText;
       if (~isempty(extraText));
          gjf_text = [gjf_text,char([10]),extraText];
       end
    end
    writeGjf( gjf_file, gjf_text );
    [origDir, tempDir] = mkScratchDir( gaussianPath ); %cd's to tempDir
    movefile( [origDir,'\',gjf_file], tempDir );
    
    setenv('GAUSS_EXEDIR', fragment.gaussianPath);

    terminated = runGaus( fragment, jobname, origDir, tempDir );
    if terminated || ~normalTermination( [tempDir,'\',jobname,'.out'] )
        bool = 0;
        cd( origDir );
        return
    end
    
    if fragment.config.opt == 1
        cd( origDir );
        fragment.opt_geom( [tempDir, '\full.out'] );
        cd( tempDir );
    end
    toZip = moveFiles( jobname, 1, 1, 1 );
    toZip = [toZip {[jobname,'.gjf']}];%, 'full_opt_config.txt'}];


%  1 NUCLEUS CALCULATIONS

if (isempty(envTarget))
   [n1,n2] = size(fragment.H1);
   natom = fragment.natom;
   fragment.H1en = zeros(n1,n2,natom);
   
   if fragment.config.calcEn == 1 && fragment.config.opt ~= 1
      tempZip = iterateAtom( fragment, origDir, tempDir );
      if isempty( tempZip )
         bool = 0;
         return
      end
      toZip = [toZip tempZip];
   end
end
%  ZIP / CLEAN UP
    
    zip(zipFileName,toZip);

    cd(origDir);
    % cleanup files
    status = rmdir(tempDir,'s');
    count1 = 0;
    while (status ~= 1) && (count1 < 1000)
       count1 = count1 + 1;
        disp('  rmdir failed. Retrying...');
        pause(0.1);
        status = rmdir(tempDir,'s');
    end
    bool = 1;
end

function toZip = iterateAtom( fragment, origDir, tempDir )
    natom = length( fragment.config.zmat.atoms );
    toZip = cell( 1, natom );
    cd( origDir );
    cd( tempDir );
    for iatom = 1:natom
        if (~fragment.config.silent)
           disp(['doing calc for atom ',num2str(iatom)]);
        end
        jobname = ['atom',num2str(iatom)];
        atom = fragment.config.zmat.atoms{iatom};
        charge = atom.z - 1;
        
        cd( origDir );
        gjf_text = buildGjf( fragment, iatom, charge, 2 );

        writeGjf( [jobname,'.gjf'], gjf_text );
        movefile( [origDir,'\',jobname,'.gjf'], tempDir );
        cd( tempDir );

        terminated = runGaus( fragment, jobname, origDir, tempDir );
        if terminated || ~normalTermination( [tempDir,'\',jobname,'.out'] )
            toZip = {};
            cd(origDir);
            return
        end
        tempZip = moveFiles( jobname, 0, 1, 0 );
        toZip{ 1, iatom } = tempZip{1};
    end
end

function terminated = runGaus(fragment, jobname, origDir, tempDir)
    %Run Gaussain with .bat file and has a timeout built in (when used)
    %Code should be backwards compatible
    %Need to do something about var named terminated

    writeGausScript( origDir );
    startTime = clock;
    timeOut = fragment.config.timeOut; % seconds
    launchBat(fragment, [fullfile(origDir, 'runGaus.bat'),' ', ...
       fullfile(tempDir,jobname), ' &'] );
    terminated = 0;
    while exist( [tempDir, '\', jobname, '.done'], 'file' ) == 0
        pause( 0.1 );
        if timeOut ~= -1 && timeCheck( startTime, timeOut )
            %I don't think this TASKKILL will work... need to
            %figure out something that works...
            system('TASKKILL /IM g09.exe /F');
            terminated = 1;
            break
        end
    end
    delete( fullfile(origDir, 'runGaus.bat'));
end

function bool = timeCheck( start, timeOut )
    %Returns 0 or 1 for if the current time (clock) is more than timeOut 
    %seconds after start
    if timeOut == -1
        bool = 0;
    else
        bool = etime( clock, start ) > timeOut;
    end
end

function [origDir, tempDir] = mkScratchDir( gaussianPath )
    tempDir = tempname([gaussianPath,'\','Scratch']);
    mkdir(tempDir);
    origDir = cd(tempDir); % save location so can move back
end

function toZip = moveFiles( jobname, moveChk, moveF32, moveOut )
    toZip = {};
    if moveChk == 1
        movefile('temp.chk', [jobname,'.chk']);
        toZip = [toZip{:} {[jobname,'.chk']}];
    end
    if moveF32 == 1
        movefile('fort.32', [jobname,'.f32']);
        toZip = [toZip{:} {[jobname,'.f32']}];
    end
    if moveOut == 1
       if (~strcmpi('full',jobname))
          movefile('full.out', [jobname,'.out']);
       end
       toZip = [toZip{:} {[jobname,'.out']}];
    end
end

function gjf = buildGjf( fragment, bq, charge, spin, extraKeywords )
    newLine = char(10);
    if nargin < 2
        bq = 0;
    end
    if nargin < 3
        charge = fragment.config.charge;
    end
    if nargin < 4
        spin = fragment.config.spin;
    end
    if nargin < 5
       extraKeywords = {};
    end

    basisSet = fragment.config.basisSet;
    method   = fragment.config.method;
    
    
    headerObj = Header( basisSet, method, fragment.config.title );
    headerObj.link0 = {'rwf=temp.rwf' 'nosave' 'chk=temp.chk'}';
    if fragment.config.opt == 1 && bq == 0
        headerObj.route = {'opt'};
    end
    headerObj.route = {headerObj.route{:} extraKeywords{:} };
    headerObj.output = {'nosymm int=noraff iop(99/6=1)' ...
        'scf=conventional' 'symm=noint'};

    headerText = headerObj.makeHeader();
    charge_mult = [num2str(charge), ' ', num2str(spin), newLine];
    zmat_body = fragment.config.zmat.build_gjf( bq );

    gjf = [headerText, charge_mult, zmat_body];
end

function writeGjf( gjf_file, gjf_text )
    fid = fopen(gjf_file,'w');
    fwrite(fid, gjf_text);
    fclose(fid);
end

function bool = normalTermination( out_loc )
    % Opens out file and searches for "Normal termination"
    
    out_text = fileread( out_loc );
    exp = 'Normal termination';
    bool = length(regexp( out_text, exp, 'once' ));
end

function writeGausScript( dir )
    newLine = char(10);
    text = ['C:\ChemSoft\G09W\g09.exe %1.gjf %1.out', newLine, ...
        'COPY /Y NUL %1.done', newLine, 'EXIT'];

    fID = fopen( fullfile(dir,'runGaus.bat'), 'w' );
    fprintf( fID, '%s', text);
    fclose( fID );
end

function launchBat(fragment, batFile)
%//LaunchBat Run a bat file with asynchronous process control
startInfo = System.Diagnostics.ProcessStartInfo('cmd.exe', sprintf('/c "%s"', batFile));
if (isempty(fragment.config.dosVisible))
startInfo.WindowStyle = System.Diagnostics.ProcessWindowStyle.Hidden;  %// if you want it invisible
end
proc = System.Diagnostics.Process.Start(startInfo);
if isempty(proc)
    error('Failed to launch process');
end
while true
    if proc.HasExited
        fprintf('\nProcess exited with status %d\n', proc.ExitCode);
        break
    end
    fprintf('.');
    pause(.1);
end
end