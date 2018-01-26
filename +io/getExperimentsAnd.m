function data = getExperimentsAnd(varargin)
% GET EXPERIMENTS AND
% fetches experiments that meet specified requirements
% The input is passed as argument pairs and an AND operater is applied to
% all matches
%
% e.g., 
% meta = io.getExperimentsAnd('Subject', 'Ellie', 'Chamber', 'V1');

meta = io.getMetaTable;

if nargin == 0 
    return
end

assert(mod(numel(varargin), 2)==0, 'Arguments must be passed as argument pairs.')

nArgs = numel(varargin)/2;
nSessions = size(meta,1);

idx = false(nSessions, nArgs);

for iArg = 1:nArgs
    
    k = (iArg-1)*2 + 1;
    arg = varargin{k};
    val = varargin{k+1};
    
    dat = meta.(arg);
    
    switch class(val)
        
        case 'char'
            
            % look for matches
            for i = 1:nSessions
               idx(i,iArg) = any(strfind(dat{i}, val));
            end
            
        case 'double'
            
            idx(:,iArg) = dat == val;
            
            
    end
        
        
    
end

idx = all(idx,2); % apply that AND

data = meta(idx,:);
