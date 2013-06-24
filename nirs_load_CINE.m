function [CINE newDir fCINE] = nirs_load_CINE(newNIRSlocation)
CINE = [];
%load TOPO.mat
try
    [newDir fil1 ext1] = fileparts(newNIRSlocation);
    fCINE = fullfile(newDir, 'CINE.mat');
    try 
        load(fCINE);
    end                
catch  exception
    disp(exception.identifier);
    disp(exception.stack(1));    
end