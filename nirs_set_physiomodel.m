function M = nirs_set_physiomodel(M)
% Physiological Model choice
M.f     = 'nirs_fx';
M.g     = 'nirs_gx';
switch M.O.PhysioModel_Choice
    case 0 %Buxton-Friston
        M.x     = zeros(4,1);
    case 1 %Zheng-Mayhew
        M.x     = zeros(5,1);    
    case 2 %Huppert1
        M.x     = zeros(7,1);           
end
M.n     = length(M.x);
%Number of inputs of direct model
M.m = 1;
%Number of outputs
st = M.IC.include_HbR + M.IC.include_HbT + M.IC.include_HbO;
if st > 1
    st = 2;
end
M.l = st;


