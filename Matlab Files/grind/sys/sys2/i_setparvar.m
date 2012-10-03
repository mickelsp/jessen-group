%sets a parameter and sets NO if iX is a state variable
function N1 = i_setparvar(iX, N0, val)
global g_grind;
N1=N0;
if ~isempty(val)
   if iX.ispar
      if isempty(iX.ndx)
         assignin('base', g_grind.pars{iX.no}, val);
      else
         multassignin('base',sprintf('%s(%d)',g_grind.pars{iX.no},iX.ndx),val);
      end;
   elseif iX.isext
      if isempty(iX.ndx)
         assignin('base', g_grind.externvars{iX.no}.name, val);
      else
         multassignin('base',sprintf('%s(%d)',g_grind.externvars{iX.no}.name,iX.ndx),val);
      end;
   else
      N1(iX.no) = val;
   end;
end;
