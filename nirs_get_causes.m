function U = nirs_get_causes(U,Stimuli)
% 'causes' or imputs U
%---------------------------------------------------------------------------
u = length(Stimuli);
U.dt = U(1).dt;

if u == 1 && length(U(1).name) == 1
    U.name = U(1).name;
    U.u    = U(1).u(33:end,1);
else
    U.name = {};
    U.u    = [];
    for  i = 1:u
        U.u             = [U.u U(Stimuli(i)).u(33:end,1)];
        U.name{end + 1} = U(Stimuli(i)).name{1};
    end
end