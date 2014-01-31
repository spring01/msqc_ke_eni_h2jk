clear a;
clear b;
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                b(i,j,k,l) = 1000*i+100*j+10*k+l;
            end
        end
    end
end
%  b = permute(b,[2 1 4 3]);
density = [1 2 3
    4 5 6
    7 8 9];

reshape(reshape(b,9,9)*reshape(density',9,1),3,3)
bk = permute(b,[1 4 2 3]);
reshape(reshape(bk,9,9)*reshape(density',9,1),3,3)

H2j = cell(3,3);
H2k = cell(3,3);
for i=1:3
    for j=1:3
        H2j{i,j} = reshape(b(i,j,:,:), 3,3);
        H2k{i,j} = reshape(b(i,:,:,j), 3,3);
    end
end
Gv = zeros(3);
Gv2 = zeros(3);
for i=1:3
    for j=1:3
        t1 = sum(sum( density'.* H2j{i,j} ));
        t2 = sum(sum( density'.* H2k{i,j} ));
        Gv(i,j) = t1;
        Gv2(i,j) = t2;
    end
end
Gv
Gv2
