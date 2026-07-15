function param = clear_basis(param)
% param = CLEAR_BASIS(param) returns a param struct with all Krylov basis
% quantities removed.

%%
fields = {'Hshort', 'Hnext', 'Vshort', 'Vnext', 'E1', 'Em'};
param = rmfield(param,fields);
end