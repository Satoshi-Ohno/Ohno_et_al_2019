%% Normalize data by the L2 norm
function     normExpData = makeNormExpData(expDataKA, optionsKA)
for f = 1 : length(expDataKA)
    tmpExpDataKA = expDataKA(f);
    tmpNorm = norm(tmpExpDataKA.level);
    tmpExpDataKA.oriL2Norm = tmpNorm;
    tmpExpDataKA.level = tmpExpDataKA.level/tmpNorm;
    tmpExpDataKA.level = double(tmpExpDataKA.level);
    normExpData(f) = tmpExpDataKA;
end

end