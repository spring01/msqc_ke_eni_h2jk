function bool = runFragment( fragment, dataPathIn, configIn )

    if (nargin < 2)
        fragment.dataPath = 'data';
        dataPath = 'data';
    else
        fragment.dataPath = dataPathIn;
        dataPath = dataPathIn;
    end
    if (nargin < 3)
        config = Config();
        fragment.config = configIn;
    else
        fragment.config = configIn;
        config = configIn;
    end

    %
    if (config.useCache)
       [found,fragment.fileprefix] = Fragment.findCalc(dataPath,config);
    else
       found = false;
    end
    if (~found)
        createTempName( fragment );
        save([fragment.fileprefix,'_cfg.mat'], 'config' );
        bool = fragment.initializeZipData( [fragment.fileprefix,'.zip'] );
        if bool == 0
            disp('initializeZipData failed')
            return
        end
    end
    fragment.loadZipData( [fragment.fileprefix,'.zip'] );
    if (~config.useCache)
       delete([ fragment.fileprefix,'.zip'] ); 
    end

    bool = 1;
end

function createTempName( fragment )
    temp1 = tempname('a'); % makes "a\uniquestring"
    uniqueStr = temp1(3:end);
    fragment.fileprefix = [fragment.dataPath,filesep,fragment.config.title, ...
        '_',uniqueStr];
end