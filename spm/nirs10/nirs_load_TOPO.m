function [TOPO newDir fTOPO] = nirs_load_TOPO(newNIRSlocation)
TOPO = [];
%load TOPO.mat
try
    [newDir fil1 ext1] = fileparts(newNIRSlocation);
    fTOPO = fullfile(newDir, 'TOPO.mat');
    try 
        load(fTOPO);
    end                
catch  exception
    disp(exception.identifier);
    disp(exception.stack(1));    
end