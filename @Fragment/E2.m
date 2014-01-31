function res = E2(obj)
% E2



Ehf = obj.Ehf;


EelecNuc = obj.Een(1);
for iatom = 2:obj.natom
   EelecNuc = EelecNuc + obj.Een(iatom);
end

res = Ehf - obj.EKE - EelecNuc - obj.Hnuc;

end