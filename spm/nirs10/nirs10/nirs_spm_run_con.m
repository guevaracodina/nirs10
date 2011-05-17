function SPM = nirs_spm_run_con(varargin)
job = varargin{1};
SPM = varargin{2};
for i = 1:length(job.consess)
    if isfield(job.consess{i},'tcon')
        name = job.consess{i}.tcon.name;        
        STAT = 'T';
        con  = job.consess{i}.tcon.convec(:)';
        sessrep = job.consess{i}.tcon.sessrep;
    elseif isfield(job.consess{i},'tconsess')
        job.consess{i}.tconsess = job.consess{i}.tconsess; % save some typing
        name = job.consess{i}.tconsess.name;
        STAT = 'T';
        if isfield(job.consess{i}.tconsess.coltype,'colconds')
            ccond = job.consess{i}.tconsess.coltype.colconds;
            con = zeros(1,size(SPM.xX.X,2)); % overall contrast
            for cs = job.consess{i}.tconsess.sessions
                for k=1:numel(ccond)
                    if SPM.xBF.order < ccond(k).colbf
                        error(['Session-based contrast %d:\n'...
                            'Basis function order (%d) in design less ' ...
                            'than specified basis function number (%d).'],...
                            i, SPM.xBF.order, ccond(k).colbf);
                    end
                    % Index into columns belonging to the specified
                    % condition
                    try
                        cind = ccond(k).colbf + ...
                            ccond(k).colmodord*SPM.xBF.order ...
                            *SPM.Sess(cs).U(ccond(k).colcond).P(ccond(k) ...
                            .colmod).i(ccond(k).colmodord+1);
                        con(SPM.Sess(cs).col(SPM.Sess(cs).Fc(ccond(k).colcond).i(cind))) ...
                            = ccond(k).conweight;
                    catch
                        error(['Session-based contrast %d:\n'...
                            'Column "Cond%d Mod%d Order%d" does not exist.'],...
                            i, ccond(k).colcond, ccond(k).colmod, ccond(k).colmodord);
                    end
                end
            end
        else % convec on extra regressors
            con = zeros(1,size(SPM.xX.X,2)); % overall contrast
            for cs = job.consess{i}.tconsess.sessions
                nC = size(SPM.Sess(cs).C.C,2);
                if nC < numel(job.consess{i}.tconsess.coltype.colreg)
                    error(['Session-based contrast %d:\n'...
                        'Contrast vector for extra regressors too long.'],...
                        i);
                end
                ccols = numel(SPM.Sess(cs).col)-(nC-1)+...
                    (0:numel(job.consess{i}.tconsess.coltype.colreg)-1);
                con(SPM.Sess(cs).col(ccols)) = job.consess{i}.tconsess.coltype.colreg;
            end
        end
        sessrep = 'none';
    else %fcon
        name = job.consess{i}.fcon.name;
        STAT = 'F';
        try
            con  = cat(1,job.consess{i}.fcon.convec{:});
        catch
            error('Error concatenating F-contrast vectors. Sizes are:\n %s\n',... 
                   num2str(cellfun('length',job.consess{i}.fcon.convec)))
        end
        sessrep = job.consess{i}.fcon.sessrep;
    end
  
    if isfield(SPM,'Sess') && ~strcmp(sessrep,'none')
        % assume identical sessions, no check!
        nsessions=numel(SPM.Sess);
        switch sessrep
            case {'repl','replsc'}
                % within-session zero padding, replication over sessions
                cons = {zeros(size(con,1),size(SPM.xX.X,2))};
                for sess=1:nsessions
                    sfirst=SPM.Sess(sess).col(1);
                    cons{1}(:,sfirst:sfirst+size(con,2)-1)=con;
                end
                if strcmp(sessrep,'replsc')
                    cons{1} = cons{1}/nsessions;
                end
                names = {sprintf('%s - All Sessions', name)};
            case 'replna',
                % within-session zero padding, new rows per session
                cons= {zeros(nsessions*size(con,1),size(SPM.xX.X,2))};
                for sess=1:nsessions
                    sfirst=SPM.Sess(sess).col(1);
                    cons{1}((sess-1)*size(con,1)+(1:size(con,1)),sfirst-1+(1:size(con,2)))=con;
                end
                names = {sprintf('%s - All Sessions', name)};
            case 'sess',
                cons = cell(1,numel(SPM.Sess));
                names = cell(1,numel(SPM.Sess));
                for k=1:numel(SPM.Sess)
                    cons{k} = [zeros(size(con,1),SPM.Sess(k).col(1)-1) con];
                    names{k} = sprintf('%s - Session %d', name, k);
                end;
            case {'both','bothsc'}
                cons = cell(1,numel(SPM.Sess));
                names = cell(1,numel(SPM.Sess));
                for k=1:numel(SPM.Sess)
                    cons{k} = [zeros(size(con,1),SPM.Sess(k).col(1)-1) con];
                    names{k} = sprintf('%s - Session %d', name, k);
                end
                if numel(SPM.Sess) > 1
                    % within-session zero padding, replication over sessions
                    cons{end+1}= zeros(size(con,1),size(SPM.xX.X,2));
                    for sess=1:nsessions
                        sfirst=SPM.Sess(sess).col(1);
                        cons{end}(:,sfirst:sfirst+size(con,2)-1)=con;
                    end
                    if strcmp(sessrep,'bothsc')
                        cons{end} = cons{end}/nsessions;
                    end
                    names{end+1} = sprintf('%s - All Sessions', name);
                end
        end
    else
        cons = {con};
        names = {name};
    end

    % Loop over created contrasts
    %-------------------------------------------------------------------
    for k=1:numel(cons)

        % Basic checking of contrast
        %-------------------------------------------------------------------
        [c,I,emsg,imsg] = spm_conman('ParseCon',cons{k},SPM.xX.xKXs,STAT);
        if ~isempty(emsg)
            disp(emsg);
            error('Error in contrast specification');
        else
            disp(imsg);
        end

        % Fill-in the contrast structure
        %-------------------------------------------------------------------
        if all(I)
            DxCon = spm_FcUtil('Set',names{k},STAT,'c',c,SPM.xX.xKXs);
        else
            DxCon = [];
        end

        %Add degrees of freedom
        DxCon.eidf = rank(c);
        
        % Append to SPM.xCon. SPM will automatically save any contrasts that
        % evaluate successfully.
        %-------------------------------------------------------------------
        if isempty(SPM.xCon)
            SPM.xCon = DxCon;
        elseif ~isempty(DxCon)
            SPM.xCon(end+1) = DxCon;
        end
        %SPM = spm_contrasts(SPM,length(SPM.xCon));
    end
end
