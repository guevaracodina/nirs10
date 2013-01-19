function TOPO = nirs_get_contrasts_group(Z,TF,TOPO)
%Organization. Either:
%1- the user has specified contrasts (TF not empty) or
%2- contrasts need to be automatically generated (TF empty), for
%each session (sessions cannot be grouped)
try
    if ~isempty(TF.consess)
        %User-specified group level contrasts
        
    else
        Sess = [];
        Cp = [];
        if Z.FFX || Z.nS == 1
            try 
                SSxCon = TOPO.xCon; % = TOPO.SSxCon; 
            catch
                disp('No SSxCon -- perhaps you are trying to run group level contrasts with only one subject!')
                disp('This will break -- make sure you have at least 2 subjects');
            end
            ns = length(SSxCon);
            
            %build up list of contrasts and of their names
            Clist = SSxCon{1};
            Nlist = {}; for t2=1:length(SSxCon{1}), Nlist = [Nlist; SSxCon{1}(t2).name]; end
            for t1=2:ns
                for t2=1:length(SSxCon{t1})
                    if ~any(strcmp(SSxCon{t1}(t2).name,Nlist))
                        %add to the list
                        try %probably never works -- 
                            Clist{end+1} = SSxCon{t1}(t2);
                            Nlist{end+1} = SSxCon{t1}(t2).name;
                        catch %always default to this?                            
                            Clist(end+1) = SSxCon{t1}(t2);
                            Nlist{end+1} = SSxCon{t1}(t2).name;
                        end
                    end
                end
            end
            nC = length(Clist);
            xCon = Clist;
            %for each contrast, build list of available sessions,
            %and position of these contrasts in SSxCon
            for c1=1:nC
                Sess{c1} = [];
                %Cp{c1} = []; %contrast position in the SSxCon list of that session
                for t1=1:ns
                    Nlist = {}; for t2=1:length(SSxCon{t1}), Nlist = [Nlist; SSxCon{t1}(t2).name]; end
                    if any(strcmp(xCon(c1).name, Nlist)),
                        Sess{c1} = [Sess{c1} t1];
                        Cp{c1,t1} = find(strcmp(xCon(c1).name, Nlist)==1);
                    end
                end
            end
            TOPO.Clist = Clist;
            TOPO.Nlist = Nlist;
        else
            xCon = TOPO.xCon;
        end
        
        TOPO.xCon = xCon;
        TOPO.Sess = Sess;
        TOPO.Cp = Cp;
        
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Could not generate group contrasts');
end