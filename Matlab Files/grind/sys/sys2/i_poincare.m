function [poincar, ivar] = i_poincare(avar, avalue, increasing)
global g_Y g_grind t;
N0 = i_initvar;
if i_settingschanged(N0)
   i_ru(g_grind.odefile, t, g_grind.ndays, N0, 1);
end;
ivar = i_varno(avar);
poincar = zeros(size(g_Y, 1), size(g_Y, 2));
if ~isempty(ivar)
   ipoin = 0;
   if increasing
      for i = 3:size(g_Y, 1) - 2
         if (g_Y(i - 1, ivar) < avalue) && (g_Y(i, ivar) >= avalue)
            ipoin = ipoin + 1;
            try
               poincar(ipoin, :) = interp1(g_Y(i - 2:i + 2, ivar), g_Y(i - 2:i + 2, :), avalue);
            catch %#ok
               poincar(ipoin, :) = g_Y(i);
            end;
         end;
      end;
   else
      for i = 3:size(g_Y, 1) - 2
         if (g_Y(i - 1, ivar) > avalue) && (g_Y(i, ivar) <= avalue)
            ipoin = ipoin + 1;
            try
               poincar(ipoin, :) = interp1(g_Y(i - 2:i + 2, ivar), g_Y(i - 2:i + 2, :), avalue);
            catch %#ok
               poincar(ipoin, :) = g_Y(i);
            end;
         end;
      end;
   end;
   poincar = poincar(1:ipoin, :);
else
   errordlg('Unknown variable');
   error('GRIND:poincare:UnknownVar','Unknown variable');
end;
