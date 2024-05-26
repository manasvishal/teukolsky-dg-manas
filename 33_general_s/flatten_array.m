function [reshaped_array] = flatten_array(array,dg_globals)


reshaped_array=reshape(array,1,[]);
reshaped_array=reshape(reshaped_array,dg_globals.Np,dg_globals.K);


end