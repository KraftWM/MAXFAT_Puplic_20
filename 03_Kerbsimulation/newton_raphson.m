function [ZVOUT, dOOUT] = newton_raphson(ZVIN, dOIN,...                    
                                         REFSIG, REFEPS,...                
                                         DEL, para, matfun,...             
                                         verfahren,...                     
                                         maxiter,tol,alpha,...             
                                         dEEPS,... 
                                         lastschritt,...
                                         varargin)
% Newton-Raphson Verfahren f�r inkrementelle Kerbn�herung
%
% INPUT:
% ZVIN        -> aktuelle Zustandsvariablen e R^(xx,1)
%                ersten beiden Tensoren m�ssen Spannungen und plastische
%                Dehungen sein
% dOIN        -> Inkrement der Energie e R^(ntens,1)
% REF...      -> Referenzzust�nde e R^(ntens,1)
% DEL         -> elastische Nachgiebigkeit e R^(ntens x ntens)
% para        -> Modell Parameter
% matfun      -> Function handle f�r Materialfunktion
% verfahren   -> Name des Verfahren, steuert Energie- und
%                Ableitungsrechnung (str)
% maxiter     -> maximal erlaubte Anzahl an Iterationen
% tol         -> Toleranz f�r abbruchbedingung
% alpha       -> Relaxationskoef.
% dEEPS       -> erstes Inkrement der Dehnungen e R^(ntens,1)
% lastschritt -> aktueller Lastschritt (f�r fehlerausgabe) e int
% varargin    -> variabler Input f�r Energie und Ableitungsrechnung
%
%
% OUTPUT:
% ZVOUT       -> neue Zustandsvariablen
% dOOUT       -> Energieinkrement mit neuen Zustandsvariablen
% _________________________________________________________________________


% Spannungszustand (aktuell nur f�r ESZ)
ntens = 3;
ndi = 2;

% Funktionen je nach verfahren
switch verfahren
    
    case 'Neuber'
        energyfun = @neuberenergy;
        ablfun = @neuberableitung;
    case 'ESED'
        energyfun = @esedenergy;
        ablfun = @esedableitung;
    case 'ModNeuber'
        energyfun = @modneuberenergy;
        ablfun = @modneuberableitung;
    case 'UniExp'
        energyfun = @uniexpenergy;
        ablfun = @uniexpableitung;
end


% Init iteration mit elastischem Dehnungsinkrement
dEPS = dEEPS;
[ZVOUT,DEP,CEP] = matfun(ntens,ndi,dEPS,ZVIN,1,para);

% Rauslesen der Zust�nde
SIG = ZVOUT(1:ntens);
EPS = DEL * SIG + ZVOUT(1+ntens:2*ntens);
dEPS = EPS - DEL * ZVIN(1:ntens) - ZVIN(ntens+1:2*ntens);
dSIG = SIG - ZVIN(1:ntens);

% berechne Fehler
dOOUT = energyfun(SIG,dSIG,EPS,dEPS,REFSIG,REFEPS,varargin{:});
err = dOIN - dOOUT;
norm2 = sum(err.*err);

% Iterationsschleife
fiter = zeros(3,maxiter+1);           % speicher zum debuggen
depsiter = zeros(3,maxiter+1);        % Speicher zum debuggen
dsigiter = zeros(3,maxiter+1);        % Speicher zum debuggen
tol_overshot = 1e-1;                  % Toleranz grenze fuer ueberschie�en
iter = 1;

while norm2 > tol
        
    % keine endlosschleifen
    iter = iter + 1;
    if iter > maxiter + 1
        msg = ['Keine Konvergenz in Newton Verfahren, ',...
               ' Anzahl Iterationen: ', num2str(iter-1),...
               ' Lastschritt: ',num2str(lastschritt),...
               ' Fehlerquadrate:', num2str(norm2)];
		warning(msg)
        break
    end
    
    % Speicher zum Debuggen
    depsiter(:,iter) = dEPS;
    dsigiter(:,iter) = dSIG;
    fiter(:,iter) = err;
    
    % Ableitung des N�herungsvrefahrens
    [G,P] = ablfun(SIG,EPS,REFSIG,REFEPS,varargin{:});
    derr = G * CEP + P;
    
    % �nderung des Dehnungsinkrements
	ddEPS = derr\err;
    
    % d�mpfen der Schrittweite
    s = 1 - alpha;         % konstante D�mpfung
    
    % neues (relaxiertes) Dehnungsinkrement
	dEPS = dEPS + s .* ddEPS;
         
    % Integration mit neuem Dehnungsinkrement
    [ZVOUT,DEP,CEP] = matfun(ntens,ndi,dEPS,ZVIN,1,para);
    
    % Rauslesen der Zust�nde
    SIG = ZVOUT(1:ntens);
    EPS = DEL * SIG + ZVOUT(1+ntens:2*ntens);
    dSIG = SIG - ZVIN(1:ntens);

    % berechne Fehler
    dOOUT = energyfun(SIG,dSIG,EPS,dEPS,REFSIG,REFEPS,varargin{:});
    err = dOIN - dOOUT;
    norm2 = sum(err.*err);

end

% Speicher zum Debuggen
depsiter(:,iter+1) = dEPS;
dsigiter(:,iter+1) = dSIG;
fiter(:,iter+1) = err;

% Plot Fehler norm
% norm2 = sum(fiter.*fiter,1);
% figure(1), plot(norm2)

% Save error to directory
% filename = ['error_alpha=',num2str(alpha),'_Inkrement_',num2str(lastschritt),'.mat'];
% save(filename,'fiter');
end % Ende Hauptfunktion

