fs.solvehf();
lambda = 0.001;
judger = 10000000;
for index=1:10000

%     judger = index.^0.8;
    newkemodmat = fs.hl.dEkestore./fs.ll.density()./fs.ll.ke();
    newkemodmat( isnan(newkemodmat) ) = 1;
    newkemodmat = (newkemodmat + newkemodmat')./2;
%     newkemodmat( newkemodmat>judger ) = 1;
%     newkemodmat( newkemodmat<1./judger ) = 1;
    fs.ll.kemodmat = (1-lambda).*fs.ll.kemodmat + lambda.*newkemodmat;
%     fs.ll.enimodmat = fs.hl.dEenistore./fs.ll.density()./fs.ll.ENI();

%     fs.llmodel.gmodmat = fs.hlmodel.decompgstore./fs.llmodel.density()./(reshape(fs.llmodel.h2jk*reshape(fs.llmodel.density(),fs.llmodel.nbasis.^2,1),fs.llmodel.nbasis,fs.llmodel.nbasis));
    fs.ll.solvehf();
    
    disp(index);
end