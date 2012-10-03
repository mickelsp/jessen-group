%%%Function for inverse beam parameter after passing through a thin
%%%lens.
function inverse2=ABCDmatrix_ThinsLens(inverse1,f)

%Calculation
inverse2= -1/f + inverse1; %new inverse beam parameter after distance z