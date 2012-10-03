%A=drawuniform(dim,minAlfa,maxAlfa)
%
%draw a matrix (A) from uniform distributions
%
function A=drawuniform(dim1,dim2,minA,maxA)
if nargin==3
   maxA=minA;
   minA=dim2;
   dim2=dim;
end;
if isempty(minA)||isempty(maxA)||isnan(minA)||isnan(maxA)
    error('GRIND:drawuniform:ArgError','Error drawuniform: range [%g,%g] is not correct\n',minA,maxA);
end;
A=rand(dim1,dim2)*(maxA-minA)+minA;
