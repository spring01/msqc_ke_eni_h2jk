
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
a = [1111 1211 1311 1121 1221 1321 1131 1231 1331
     2111 2211 2311 2121 2221 2321 2131 2231 2331
     3111 3211 3311 3121 3221 3321 3131 3231 3331
     1112 1212 1312 1122 1222 1322 1132 1232 1332
     2112 2212 2312 2122 2222 2322 2132 2232 2332
     3112 3212 3312 3122 3222 3322 3132 3232 3332
     1113 1213 1313 1123 1223 1323 1133 1233 1333
     2113 2213 2313 2123 2223 2323 2133 2233 2333
     3113 3213 3313 3123 3223 3323 3133 3233 3333];
 
reshape(a,3,3,3,3) - permute(b,[1 4 2 3])