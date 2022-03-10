function [SIG, EPS, EPSP, ALPHAI, R, P, ZVAR1, EZVAR1] = ...
          pseudostresslang_ZV_auslesen(ZVAR,material,para,epara,M,eM,numink)
% Funktion zum auslesen der Zustandsvariablen aus dem Pseudo Stress 
% approach nach Lang      
%
% INPUT:
% ZVAR              - Zustandsvariablen
% material          - name des Plastizitaetsmodells
% para              - Parameter des Plastizit�smodells !!! an die para
%                     k�nnten die epara geh�ngt sein
% epara             - Parameter des Strukturmodells
% M                 - Anzahl Backstresstensoren
%
% OUTPUT:
% Spannungen, Dehnungen, plastische Dehnungen ....
%
%__________________________________________________________________________

% !!!!!!!!!!!!!!!! nur ESZ
ntens = 3;
ndi = 2;


% =========================================================================S
% Herrauslesen der lokalen Gr��en je nach Material
switch material
    % ---------------------------------------------------------------------
    case 'Chaboche'
        
        % =================================================================
        % Herrauslesen gr��en
        % Abbildungen
        P = set_maps(ntens,ndi);
        % Steifigkeitsmatrix
        CEL = elast_steifigkeit(para(1),para(2),ntens,ndi);
        DEL = elast_nachgiebigkeit(para(1),para(2),ntens,ndi);
        % pseudo gesamtdehnung 
        EEPS =  ZVAR(1:ntens,:);
        % plastische Dehnung
        EPSP = ZVAR(ntens+1:2*ntens,:);
        % Pseudo (elastische9 Spannungen
        ESIG = CEL * (EEPS - EPSP);
        % Pseudo Spannungs Deviator
        ES = P * ESIG;
        % pseudo backstress
        EALPHAI = reshape(ZVAR(2*ntens+1:(eM+2)*ntens,:),ntens,eM,numink);        
        % "realer" Backstress
        ALPHAI = reshape(ZVAR((eM+2)*ntens+3:(eM+M+2)*ntens+2,:),ntens,M,numink);
%         ALPHAI = reshape(ZVAR((M+2)*ntens+2:(2*M+2)*ntens+1,:),ntens,M,numink);
        % ... Radius nicht in ZVAR
%         % plastische Bogenl�nge
%         P = ZVAR((M+2)*ntens+1,:);
%         % Ausgeben der Ged�chtnissfl�che anstatt R
%         q = para(3);
%         b = para(4); if b==0, b=1e-40; end 
%         r0 = para(5+2*M);
%         Qinf = q/b;
%         eq = epara(3);
%         eb = epara(4); if eb==0, eb=1e-40; end 
%         er0 = epara(5+2*M);
%         eQinf = eq/eb;
%         R = r0 + Qinf*(1-exp(-b.*P));
%         ER = er0 + eQinf*(1-exp(-eb.*P));
        % ... Radius in ZVAR
        % plastische Bogenl�nge
        P = ZVAR((eM+2)*ntens+2,:);
        % Radien
        R = ZVAR((eM+M+2)*ntens+3,:);
        ER = ZVAR((eM+2)*ntens+1,:);
        % =================================================================
        % berechnen des "realen" Spannungsdeviators und so

        % Backstresstensoren
        A = reshape(sum(ALPHAI,2),ntens,numink);
        EA = reshape(sum(EALPHAI,2),ntens,numink);

        % Spannungsdeviator
        S = A + R./ER .* (ES - EA);

        % Spannungen
        SIG = [2,1,0;1,2,0;0,0,1]* S;
        
        % Dehungen
        EPS = EPSP + DEL * SIG;
        % =================================================================
        % Endwerte
        
        % Strukturmodell
        EZVAR1 = [EEPS(:,end);EPSP(:,end);reshape(EALPHAI(:,:,end),3*eM,1);ER(end);P(end)];
        % Materialmodell
        ZVAR1 = [SIG(:,end);EPSP(:,end);reshape(ALPHAI(:,:,end),3*M,1);R(end);P(end)];
    
    % ---------------------------------------------------------------------    
    case 'OhnoWang'
        
        % =================================================================
        % Herrauslesen gr��en
        % Abbildungen
        P = set_maps(ntens,ndi);
        % Steifigkeitsmatrix
        CEL = elast_steifigkeit(para(1),para(2),ntens,ndi);
        DEL = elast_nachgiebigkeit(para(1),para(2),ntens,ndi);
        % pseudo gesamtdehnung 
        EEPS =  ZVAR(1:ntens,:);
        % plastische Dehnung
        EPSP = ZVAR(ntens+1:2*ntens,:);
        % Pseudo (elastische9 Spannungen
        ESIG = CEL * (EEPS - EPSP);
        % Pseudo Spannungs Deviator
        ES = P * ESIG;
        % pseudo backstress
        EALPHAI = reshape(ZVAR(2*ntens+1:(eM+2)*ntens,:),ntens,eM,numink);
        % plastische Bogenl�nge
        P = ZVAR((eM+2)*ntens+1,:);
        % "realer" Backstress
        ALPHAI = reshape(ZVAR((eM+2)*ntens+2:(eM+M+2)*ntens+1,:),ntens,M,numink);
        % Ausgeben der Ged�chtnissfl�che anstatt R
        R = para(3+3*M) * ones(1,numink);
        ER = epara(3+3*eM) * ones(1,numink);
        % =================================================================
        % berechnen des "realen" Spannungsdeviators und so

        % Backstresstensoren
        A = reshape(sum(ALPHAI,2),ntens,numink);
        EA = reshape(sum(EALPHAI,2),ntens,numink);

        % Spannungsdeviator
        S = A + R./ER .* (ES - EA);

        % Spannungen
        SIG = [2,1,0;1,2,0;0,0,1]* S;
        
        % Dehungen
        EPS = EPSP + DEL * SIG;
        % =================================================================
        % Endwerte
        
        % Strukturmodell
        EZVAR1 = ZVAR(1:3*(2+eM)+1,end);
        % Materialmodell
        ZVAR1 = [SIG(:,end);EPSP(:,end);reshape(ALPHAI(:,:,end),3*M,1);P(end)];
        
        
        
    % ---------------------------------------------------------------------    
    case 'KarimOhno'
        
        % =================================================================
        % Herrauslesen gr��en
        % Abbildungen
        P = set_maps(ntens,ndi);
        % Steifigkeitsmatrix
        CEL = elast_steifigkeit(para(1),para(2),ntens,ndi);
        DEL = elast_nachgiebigkeit(para(1),para(2),ntens,ndi);
        % pseudo gesamtdehnung 
        EEPS =  ZVAR(1:ntens,:);
        % plastische Dehnung
        EPSP = ZVAR(ntens+1:2*ntens,:);
        % Pseudo (elastische9 Spannungen
        ESIG = CEL * (EEPS - EPSP);
        % Pseudo Spannungs Deviator
        ES = P * ESIG;
        % pseudo backstress
        EALPHAI = reshape(ZVAR(2*ntens+1:(eM+2)*ntens,:),ntens,eM,numink);
        % plastische Bogenl�nge
        P = ZVAR((eM+2)*ntens+1,:);
        % "realer" Backstress
        ALPHAI = reshape(ZVAR((eM+2)*ntens+2:(eM+M+2)*ntens+1,:),ntens,M,numink);
        % Ausgeben der Ged�chtnissfl�che anstatt R
        R = para(3+3*M) * ones(1,numink);
        ER = epara(3+3*eM) * ones(1,numink);
        % =================================================================
        % berechnen des "realen" Spannungsdeviators und so

        % Backstresstensoren
        A = reshape(sum(ALPHAI,2),ntens,numink);
        EA = reshape(sum(EALPHAI,2),ntens,numink);

        % Spannungsdeviator
        S = A + R./ER .* (ES - EA);

        % Spannungen
        SIG = [2,1,0;1,2,0;0,0,1]* S;
        
        % Dehungen
        EPS = EPSP + DEL * SIG;
        % =================================================================
        % Endwerte
        
        % Strukturmodell
        EZVAR1 = ZVAR(1:3*(2+eM)+1,end);
        % Materialmodell
        ZVAR1 = [SIG(:,end);EPSP(:,end);reshape(ALPHAI(:,:,end),3*M,1);P(end)];
        
    % ---------------------------------------------------------------------    
    case 'Jiang'
        
        % =================================================================
        % Herrauslesen gr��en
        % Abbildungen
        P = set_maps(ntens,ndi);
        % Steifigkeitsmatrix
        CEL = elast_steifigkeit(para(1),para(2),ntens,ndi);
        DEL = elast_nachgiebigkeit(para(1),para(2),ntens,ndi);
        % pseudo gesamtdehnung 
        EEPS =  ZVAR(1:ntens,:);
        % plastische Dehnung
        EPSP = ZVAR(ntens+1:2*ntens,:);
        % Pseudo (elastische9 Spannungen
        ESIG = CEL * (EEPS - EPSP);
        % Pseudo Spannungs Deviator
        ES = P * ESIG;
        % pseudo backstress
        EALPHAI = reshape(ZVAR(2*ntens+1:(eM+2)*ntens,:),ntens,eM,numink);
        % plastische Bogenl�nge
        P = ZVAR((eM+2)*ntens+1,:);
        % pseudo ged�chnissfl�che
        ERM = ZVAR((eM+2)*ntens+2,:);
        % "realer" Backstress
        ALPHAI = reshape(ZVAR((eM+2)*ntens+3:(eM+M+2)*ntens+2,:),ntens,M,numink);
        % "reale" ged�chtnissfl�che
        RM = ZVAR((eM+M+2)*ntens+3,:);
        % Ausgeben der Ged�chtnissfl�che anstatt R
        R = RM;

        % =================================================================
        % berechnen des "realen" Spannungsdeviators und so

        % Radien Flie�fl�chen
        ak = para(5);
        ck = para(6);
        k1 = para(7);
        eak = epara(5);
        eck = epara(6);
        ek1 = epara(7);
        K = k1 .* (ones(size(RM)) + ak .* exp( ck .* RM ) );
        EK = ek1 .* (ones(size(ERM)) + eak .* exp( eck .* ERM ) );

        % Backstresstensoren
        A = reshape(sum(ALPHAI,2),ntens,numink);
        EA = reshape(sum(EALPHAI,2),ntens,numink);

        % Spannungsdeviator
        S = A + K./EK .* (ES - EA);

        % Spannungen
        SIG = [2,1,0;1,2,0;0,0,1]* S;
        
        % Dehungen
        EPS = EPSP + DEL * SIG;
        
        % =================================================================
        % Endwerte
        
        % Strukturmodell
        EZVAR1 = ZVAR(1:3*(2+eM)+2,end);
        % Materialmodell
        ZVAR1 = [SIG(:,end);EPSP(:,end);reshape(ALPHAI(:,:,end),3*M,1);P(end);RM(end)];
        
    case 'OWT'
        
        % =================================================================
        % Herrauslesen gr��en 
        % Abbildungen
        P = set_maps(ntens,ndi);
        % Steifigkeitsmatrix
        CEL = elast_steifigkeit(para(1),para(2),ntens,ndi);
        DEL = elast_nachgiebigkeit(para(1),para(2),ntens,ndi);
        % pseudo gesamtdehnung 
        EEPS =  ZVAR(1:ntens,:);
        % plastische Dehnung
        EPSP = ZVAR(ntens+1:2*ntens,:);
        % Pseudo (elastische9 Spannungen
        ESIG = CEL * (EEPS - EPSP);
        % Pseudo Spannungs Deviator
        ES = P * ESIG;
        % pseudo backstress
        EALPHAI = reshape(ZVAR(2*ntens+1:(eM+2)*ntens,:),ntens,eM,numink);
        % plastische Bogenl�nge
        P = ZVAR((eM+2)*ntens+1,:);
        % zus�tzlicher Radius FF Strukturmodell
        EQ = ZVAR((eM+2)*ntens+2,:);
        % !!! Zus�tzliche Parameter nicht mit Auslesen (werden f�r die
        % Auswertung nicht gebraucht
        % "realer" Backstress
        ALPHAI = reshape(ZVAR((3+eM+(ntens+1)/2)*ntens+5 : (3+eM+M+(ntens+1)/2)*ntens+4,: ),ntens,M,numink);
        % zus�tzlicher Radius FF Materialmodell
        Q = ZVAR((3+eM+M+(ntens+1)/2)*ntens+5,:);
        % Radien der Flie�fl�chen Struktur- & Materialmodell
        R = para(11+3*M) * ones(1,numink) + Q;
        ER = epara(11+3*eM) * ones(1,numink) + EQ;
        
        % =================================================================
        % berechnen des "realen" Spannungsdeviators und so

        % Backstresstensoren
        A = reshape(sum(ALPHAI,2),ntens,numink);
        EA = reshape(sum(EALPHAI,2),ntens,numink);

        % Spannungsdeviator
        S = A + R./ER .* (ES - EA);

        % Spannungen
        SIG = [2,1,0;1,2,0;0,0,1]* S;
        
        % Dehungen
        EPS = EPSP + DEL * SIG;
        % =================================================================
        % Endwerte
        
        % Strukturmodell
        EZVAR1 = ZVAR(1:(5+eM)*3+4,end);
        % Materialmodell
        ZVAR1 = [SIG(:,end);EPSP(:,end);reshape(ALPHAI(:,:,end),3*M,1);P(end);ZVAR((2+M)*ntens+2:(3+M+(ntens+1)/2)*ntens+4,end)];
        
     otherwise
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
end
      
end % Ende Funktion