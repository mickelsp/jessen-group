%gets a parameter and gets [] if iX is a state variable
function val=i_getparvar(iX)
global g_grind;
if iX.ispar
   if isempty(iX.ndx)
      val = evalin('base', g_grind.pars{iX.no});
   else
      val =evalin('base',sprintf('%s(%d)',g_grind.pars{iX.no},iX.ndx));
   end;
elseif iX.isext
   if isempty(iX.ndx)
      val = evalin('base', g_grind.externvars{iX.no}.name);
   else
      val =evalin('base',sprintf('%s(%d)', g_grind.externvars{iX.no}.name,iX.ndx));
   end;
else
   val=[];
end;
