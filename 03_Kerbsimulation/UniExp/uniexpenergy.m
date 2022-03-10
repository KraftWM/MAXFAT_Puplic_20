function dTSED = uniexpenergy(SIG,dSIG,EPS,dEPS,REFSIG,REFEPS,Cq)
% Was macht Funktion?
    dTSED = (1+Cq) * (SIG - REFSIG) .* dEPS + ... 
            (1-Cq) * (EPS - REFEPS) .* dSIG;
end