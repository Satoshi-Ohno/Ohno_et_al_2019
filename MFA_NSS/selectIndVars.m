%%     Select independent variables
function [idIndVars, idRedundantConstr,  Cx, C0, idInvalidIndVars] = ...
    selectIndVars(A, isInputFullIndVarsConstr, idInputIndVars, idInputRedundantConstr)

if nargin<=3
    idInputRedundantConstr= [];
end
if nargin<=2
    idInputIndVars = [];
end
if nargin <= 1
    isInputFullIndVarsConstr = false;
end

%% Identify redundanct variables
Aori = A;
isVarsOri = any(Aori,1);
idVarsOri = find(isVarsOri);
nVarsOri = length(idVarsOri);
idNotVarsOri = find(~isVarsOri);
if isInputFullIndVarsConstr
    inputIndVarsOri = idInputIndVars;
    isInputIndVarsOri = false(1,size(Aori,2));
    isInputIndVarsOri(inputIndVarsOri) = true;
    isInputIndVars = isInputIndVarsOri(isVarsOri);
    idInputIndVars = find(isInputIndVars);
end

%% Remove redundanct variables
A = Aori(:,isVarsOri);
[nConstr, nVars] = size(A);
nInputIndVars = length(idInputIndVars);

%% Add constraints
if ~isInputFullIndVarsConstr
    if  nInputIndVars>=1
        A = [A; zeros(length(idInputIndVars), nVars)];
        A(end-nInputIndVars +1: end, idInputIndVars) = ...
            A(end-nInputIndVars +1: end, idInputIndVars) + eye(length(idInputIndVars));
    end
end

%% Identify indepndent constraints
if isInputFullIndVarsConstr
    isRedundantConstr = false(size(A,1),1);
    isRedundantConstr(idInputRedundantConstr) = true;
    idRedundantConstr = find(isRedundantConstr);
    idDepConstr = find(~isRedundantConstr);
else
    [idRedundantConstr, idDepConstr]=identifyIndVar(A');
    isRedundantConstr = false(size(A,1),1);
    isRedundantConstr(idRedundantConstr) = true;
end

%% Check wheter input parameters are independet or not
if isInputFullIndVarsConstr
else
    id = find(idRedundantConstr>nConstr);
    idInvalidIndVars = zeros(length(id),1);
    for i = 1 :length(id)
        idInvalidIndVars(i) = idInputIndVars(idRedundantConstr(id(i))-nConstr);
    end
    
    loc = ismember(idInputIndVars, idInvalidIndVars);
    idInputIndVars = idInputIndVars(~loc);
    nInputIndVars = length(idInputIndVars);
    if  ~isempty(idInputIndVars)
        A = A(1:nConstr, :);
        A = [A; zeros(length(idInputIndVars), nVars)];
        A(end-nInputIndVars +1: end, idInputIndVars) = ...
            A(end-nInputIndVars +1: end, idInputIndVars) + eye(length(idInputIndVars));
    end
    
    [idRedundantConstr, idDepConstr]=identifyIndVar(A');
    [idRedundantConstrTest, idDepConstrTest]=identifyIndVar(A(idDepConstr,:)');
    isRedundantConstr = false(size(A,1),1);
    isRedundantConstr(idRedundantConstr) = true;
end
isRedundantConstrOri = isRedundantConstr(1:nConstr);
idDepConstrOri = idDepConstr(idDepConstr<=nConstr);

%% Identify indepdent parameters
if isInputFullIndVarsConstr
    idNewIndVars = idInputIndVars;
    idIndVars = idInputIndVars;
    idDepVars = find(~isInputIndVars);
else
    [idNewIndVars, idDepVars]=identifyIndVar(A);
    if nInputIndVars == 0
        idIndVars = idNewIndVars;
    else
        idIndVars = [idInputIndVars; idNewIndVars];
        idDepVars = setdiff(idDepVars, idInputIndVars);
    end
end
idIndVars = sort(idIndVars); 
nIndVars = length(idIndVars);

%% Make matrix to calculate all parameters from idependent parameters
% Ax = b;
% Adep * xdep + Aind * xind = b;
% xdep = Adep^-1 * (b-Aind*xind) = Adep^-1*b  - Adep^-1*Aind*xind;

Cdep = inv(A(idDepConstrOri,idDepVars));

Cx = sparse(nVars,length(idIndVars));
Cx(1:size(Cdep,1), :) = -Cdep*A(idDepConstrOri,idIndVars);
Cx(size(Cdep,1)+1:end, :) = eye(length(idIndVars));

C0 = sparse(nVars, length(idDepConstrOri));
C0(1:size(Cdep,1), :) = Cdep;
C0(size(Cdep,1)+1:end, :) = 0;

Cx= round(Cx*10^12)/10^12;
C0= round(C0*10^12)/10^12;

[~,idSort]=sort([colVec(idDepVars); colVec(idIndVars)]);
Cx = Cx(idSort,:);
C0 = C0(idSort,:);

%% Modify
% redundant constrains and variables are considerd

isIndVarsOri = false(1,nVarsOri);
isIndVarsOri(idVarsOri(idIndVars)) = true;
idIndVarsOri = find(isIndVarsOri);


Cxori = sparse(numel(isVarsOri), nnz(~isVarsOri)+nIndVars);
Cxori(~isVarsOri, 1:nnz(~isVarsOri)) = ...
    Cxori(~isVarsOri, 1:nnz(~isVarsOri)) + speye(nnz(~isVarsOri));
Cxori(isVarsOri, nnz(~isVarsOri)+(1:nIndVars)) = ...
    Cxori(isVarsOri, nnz(~isVarsOri)+(1:nIndVars)) + Cx;

C0ori = sparse(numel(isVarsOri), numel(~isRedundantConstrOri));
C0ori(isVarsOri, ~isRedundantConstrOri) = C0ori(isVarsOri, ~isRedundantConstrOri) + C0;

[idIndVarsOri, idSort] = sort([idNotVarsOri, idIndVarsOri]);
Cxori = Cxori(:,idSort);

Cx = Cxori;
C0 = C0ori;
idIndVars = colVec(idIndVarsOri);


end

%% Identify independent variabls
function  [idIndVar, idDepVar]=identifyIndVar(A)

rrefA = rref(A, 10^-9);  

tmpColId = zeros(size(rrefA,1),1);
for r = 1 : size(rrefA,1)
    idCol =  find(rrefA(r,:),1);
    if isempty(idCol)
        continue
    end
    tmpColId(r) = idCol;
end
idIndRow = find(tmpColId~=0);
tmpColId = [tmpColId(idIndRow);size(A,2)+1];
diffCol = diff(tmpColId);

idIndVar = [];
for r = 1 : length(tmpColId)-1
    if diffCol(r) == 1
        continue
    end
    idIndVar= [idIndVar, tmpColId(r)+(1:diffCol(r)-1)];
end
idIndVar = colVec(idIndVar);
idDepVar = colVec(setdiff(1:size(A,2), idIndVar));

end

 