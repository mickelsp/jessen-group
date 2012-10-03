%%%Function for inverse beam parameter after some longitudinal
%%%distance from where the original inverse beam parameter is defined
function inverse2=ABCDmatrix_FreeSpace(inverse1,z)

%Computation
inverse2=inverse1./(1+z.*inverse1); %new inverse beam parameter after distance z