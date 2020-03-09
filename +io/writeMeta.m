function writeMeta(new_meta, overwrite)

if nargin < 2
    overwrite = 1;
end

[old_meta, meta_file] = io.getMetaTable();
meta = old_meta; % copy to be written

nNew = size(new_meta,1);

for iSession = 1:nNew
    
    fields = meta.Properties.VariableNames;
    
    % check if the session already exists
    old_ix = find(strcmp(old_meta.Directory, new_meta.Directory{iSession}));
    
    if isempty(old_ix)
        
        fprintf('Session %d has no matches. Appending to meta\n', iSession)
        
        meta = [meta; new_meta(iSession,:)];
        
    else
        
        old_ = meta(old_ix,:);
        new_ = new_meta(iSession,:);
        
        nVar = numel(fields);
        mismatch = false(1,nVar);
        
        for iVar = 1:nVar
            
            if ~isa(old_.(fields{iVar}), class(new_.(fields{iVar})))
                mismatch(iVar) = true;
            else
                switch class(old_.(fields{iVar}))
                    
                    case 'cell'
                        
                        if iscell(new_.(fields{iVar}))
                            
                            if (numel(old_.(fields{iVar}){1}) ~= numel(new_.(fields{iVar}){1})) || ~all(old_.(fields{iVar}){1} == new_.(fields{iVar}){1})
                                mismatch(iVar) = true;
                            end
                            
                            if mismatch(iVar) % if mismatch flagged, check for nans
                                if (numel(old_.(fields{iVar}){1})==1) && isnan(old_.(fields{iVar}){1}) && (numel(new_.(fields{iVar}){1})==1) && isnan(new_.(fields{iVar}){1})
                                    mismatch(iVar) = false;
                                end
                            end
                        else
                            new_.(fields{iVar}) = {new_.(fields{iVar})};
                        end
                    case 'double'
                        
                        if ~all(old_.(fields{iVar}) == new_.(fields{iVar}))
                            mismatch(iVar) = true;
                        end
                        
                        if mismatch(iVar) % if mismatch flagged, check for nans
                            if isnan(old_.(fields{iVar})) && isnan(new_.(fields{iVar}))
                                mismatch(iVar) = false;
                            end
                        end
                end
            end
            
            
            
        end
        
        if any(mismatch)
            
            fprintf('Session %d exists and does not match for %d variables\n', iSession, sum(mismatch))
            for i = find(mismatch(:)')
                fprintf('%s\n', fields{i})
            end
            
            disp('Old Entry:')
            disp(old_)
            
            disp('New Entry:')
            disp(new_)
            
            switch overwrite
                
                case 0
                    disp('Overwrite is false. Skipping.')
                    
                case 1
                    
                    reply = input('Do you want to overwrite? (1 or 0)');
                    
                    if reply
                        meta(old_ix,:) = new_;
                    else
                        disp('skipping')
                    end
                    
                case 2
                    meta(old_ix,:) = new_;
            end
        else
            disp('No changes to this session. Skipping')
        end
        
    end
    
end

if overwrite==2
    reply = true;
else
    reply = input('Do you want to overwrite the meta data file? (1 or 0)');
end

if reply
    C = table2cell(meta); C = [meta.Properties.VariableNames; C];
    
    disp('Saving.')
    xlswrite(meta_file, C);
    disp('Meta File Updated!')
end





