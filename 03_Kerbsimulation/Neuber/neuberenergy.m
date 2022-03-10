function dTSED = neuberenergy(SIG,dSIG,EPS,dEPS,REFSIG,REFEPS)
% Was macht Funktion
    dTSED = (SIG - REFSIG) .* dEPS + ... 
            (EPS - REFEPS) .* dSIG;
end