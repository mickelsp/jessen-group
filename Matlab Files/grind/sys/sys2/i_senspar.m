function i_senspar(flag, par1, par2)
global g_grind g_senspar g_Y t;
switch flag
 case '-init'
   evalin('base','global g_senspar;');
   g_senspar = [];
   g_senspar.n=par2;
   g_senspar.pars = {};
   for i = 1:length(g_grind.pars);
      par.name = g_grind.pars{i};
      par.val = evalin('base', g_grind.pars{i});
      par.minval = par.val * (1 - par1);
      par.maxval = par.val * (1 + par1);
      par.vals = zeros(g_senspar.n,1) * NaN;
      g_senspar.pars{i} = par;
   end;
 case '-montecarlo'
   oldstep = g_grind.tstep;
   i_parcheck;
   if isnan(oldstep)
      g_grind.tstep = 10;
   end;
   g_senspar.Ys=zeros(g_senspar.n,g_grind.statevars.dim*g_grind.tstep)*NaN;
   g_grind.solver.opt.OutputFcn = [];
   h = waitbar(0, 'Simulating');
   for j = 1:g_senspar.n
      waitbar(j / g_senspar.n, h);
      set(h,'name',sprintf('Run %d/%d',j,g_senspar.n));
      for i = 1:length(g_senspar.pars)
         par = g_senspar.pars{i};
         v = par.minval + rand(1) * (par.maxval - par.minval);
         g_senspar.pars{i}.vals(j) = v;
         assignin('base',par.name, v);
      end;
      N0 = i_initvar;
      i_ru(g_grind.odefile, t, g_grind.ndays, N0, 1);
      y = g_Y(2:end, :)';
      g_senspar.Ys(j, :) = y(:)';
   end;
   close(h);
   for i = 1:length(g_senspar.pars)
      par = g_senspar.pars{i};
      assignin('base', par.name,par.val);
   end;
   g_grind.tstep = oldstep;
case '-matrix'
   for p=1:length(g_senspar.pars)
      for k=1:size(g_senspar.Ys,2)
         v=polyfit(g_senspar.pars{p}.vals,g_senspar.Ys(:,k),1);
         g_senspar.sensmatrix(p,k)=v(1)/mean(g_senspar.pars{p}.vals)*...
            (g_senspar.pars{p}.maxval-g_senspar.pars{p}.minval);
      end;
   end;
otherwise
   fprintf('i_senspar: "%s" not implemented\n',flag);
end;


