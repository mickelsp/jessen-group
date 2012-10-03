function v=initval(astatevar,par)
if nargin==0
    val;
elseif (nargin==1) && (nargout==1)  
   v=evalin('base',astatevar);
elseif (nargin==1)
   val(astatevar);
elseif (nargin==2)
   val(astatevar,par);
end;