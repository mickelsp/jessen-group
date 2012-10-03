function g_Jacobian = i_calcjac(donumerical,niters, g_l_N0)
global g_grind g_t; %#ok g_t is needed  
g_l=g_grind.statevars.dim;
g_Jacobian = zeros(g_l);
if ~donumerical
   if ~isempty(g_grind.pars)
      eval(i_globalstr(g_grind.pars));
   end;
   for g_l_i = 1:g_l
      eval(sprintf('%s=g_l_N0(%d);',i_statevars_names(g_l_i),g_l_i));
   end;
   if ~isempty(g_grind.externvars)
     for g_l_m = 1:length(g_grind.externvars)
        eval(sprintf('%s=externvar(g_grind.externvars{%d}.data, %s,g_t,g_grind.externvars{%d}.options);',g_grind.externvars{g_l_m}.name, ...
           g_l_m, g_grind.externvars{g_l_m}.default, g_l_m));
      end;
   end;
   for g_l_i = 1:g_l
      for g_l_j = 1:g_l
         g_Jacobian(g_l_i,g_l_j)=eval(g_grind.Jacobian{g_l_i,g_l_j});
      end
   end
   if g_grind.solver.backwards
      g_Jacobian = -g_Jacobian;
   end;
else
   if g_grind.solver.backwards
      if g_grind.solver.isdiffer
         odefile = 'i_backdiffer';
     else
         odefile = 'i_back';
      end;
   else
      odefile = g_grind.odefile;
   end;
   Nres1 = g_l_N0;
   for k = 1:niters
      Nres1 = feval(str2func(odefile), 1, Nres1);
   end;
   delta = 0.001;
   for i = 1:g_l
      N0_1 = g_l_N0;
      N0_1(i) = g_l_N0(i) + delta;
      for k = 1:niters
         N0_1 = feval(odefile, 1, N0_1);
      end;
      Nres2 = N0_1;
   %   for j = 1:g_l
         g_Jacobian(:, i) = (Nres2 - Nres1) / delta;
   %   end
   end
end
