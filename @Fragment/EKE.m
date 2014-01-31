function res = EKE(obj)
% Kinetic energy

res = sum(sum( obj.density.*obj.KE ) );
end

