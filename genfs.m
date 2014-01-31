% clear classes;
clear;
filebasename = 'met';
fullfilename = strcat('C:\working\matlab\msqc_fith1g\fitting\', filebasename, '.mat');

fs = FitSingle(fullfilename,7);

% iniparvec = ones(1,fs.numpars);
% iniparvec = x + 0.01.*(rand(1,fs.numpars) - 0.5);
% lowlimits = [];
% highlimits = [];
% options = optimset('MaxFunEvals',1e6,'Display','iter','TolX',1e-16,'TolFun',1e-16,'algorithm','levenberg-marquardt','MaxIter',1e4);
% options = optimset('MaxFunEvals',1e6,'Display','iter','TolX',1e-14,'TolFun',1e-14,'MaxIter',1e4);
% options = optimset('MaxFunEvals',1e6,'Display','iter','algorithm','levenberg-marquardt','MaxIter',1e4);
% x = lsqnonlin(@fs.err, iniparvec, lowlimits, highlimits, options);
% x = fsolve(@fs.err,iniparvec,options);



% fs = cell(1,10);
% for i=1:10
%     fs{i} = FitSingle(fullfilename,i);
%     iniparvec = ones(1,fs{i}.numpars);
%     fstmp = fs{i};
%     x(i,1:fstmp.numpars) = lsqnonlin(@fstmp.err, iniparvec, lowlimits, highlimits, options);
% end


