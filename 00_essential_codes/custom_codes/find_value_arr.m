function [row_idx,col_idx]=find_value_arr(arr,value)

%% MV: this code finds the given value inside an array. Since for the dG grid
%% we will have two indices, I choose to use the first index i.e., the outflow
%% point in the previous grid (Np,col_idx). I decided to do this because we found out
%% that the outflow points has faster convergence.


[Np,K]=size(arr);
location_interface=1;
value_arr=zeros(size(arr))+value;

tol=min(min(abs(value_arr-arr)));

[row_idx,col_idx]=find( abs((arr-value))<=tol );


if location_interface==1 
    row_idx=Np;
    col_idx=col_idx(1);

end

col_idx=col_idx(1);

end