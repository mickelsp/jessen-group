function i_fillpars(plusstatvars)
global g_grind;
if ~isempty(g_grind.model)
   if nargin == 0
      plusstatvars = 0;
   end;
   pars = {};
   gmodel = cell(1, size(g_grind.model, 2));
   ig = 1;
   infun = 0;
   %strip matlabfunctions out of g_grind.model
   for i = 1:size(g_grind.model, 2);
      s = g_grind.model{i};
      if strncmp(s, 'function ', 9)
         infun = 1;
      end;
      if ~infun && ~strncmp(s, 'defextern ',10)&&~strncmp(s, 'defpermanent ',13)
         gmodel{ig} = s;
         ig = ig + 1;
      end;
      if strncmp(s, 'return', 6)
         infun = 0;
      end;
   end;
   gmodel = gmodel(1:ig - 1);
   for i = 1:ig - 1
      f = strfind(gmodel{i}, ' (');
      for j = length(f):-1:1
         k = 1;
         while (f(j) - k > 1) && (gmodel{i}(f(j) - k) == ' ')
            k = k + 1;
         end;
         if (f(j) > 1) && isempty(strfind('*/+-^', gmodel{i}(f(j)-k)))
            gmodel{i} = [gmodel{i}(1:f(j) - k) gmodel{i}(f(j) + 1:length(gmodel{i}))];
         end;
      end;
   end;
   ifilter={'pi','inf','Inf','nan','NaN','else','if','end','while','for',...
      'switch','otherwise','case','global','function','return','definepars','defextern','defpermanent'};
   gstat = cell(1, size(gmodel, 2));
   ip1 = 1;
   if plusstatvars
      for i = 1:size(gmodel, 2)
         if ~isempty(strtrim(gmodel{i})) && (gmodel{i}(1) ~= '%')
            pars = [pars ; i_symvar(gmodel{i}, ifilter)];
            i2=strfind(gmodel{i}, '=');
            if ~isempty(i2)
               j1=strfind(gmodel{i}, '''');
               if ~isempty(j1) && (j1(1) > i2(1))
                  j1 = [];
               end;
               j2 = strfind(gmodel{i}, '(');
               if ~isempty(j2) && (j2(1) > i2(1))
                  j2 = [];
               end;
               if (isempty(j1) && isempty(j2))
                  gstat{ip1} = i_symvar(gmodel{i}(1:i2), ifilter);
                  ip1 = ip1 + 1;
               end;
            end;
         end;
      end;
   else
      for i = 1:size(gmodel, 2)
         if ~isempty(strtrim(gmodel{i}))
            if(gmodel{i}(1) ~= '%')
            pars = [pars ; i_symvar(strrep(gmodel{i},'rednoise(t,','rednoise('), ifilter)];
            i2=strfind(gmodel{i}, '=');
            if ~isempty(i2)
               i3=strfind(gmodel{i}(1:i2),'''');
               if ~isempty(i3)
                  i2 = i3;
               end;
               i3 = strfind(gmodel{i}(1:i2), '(');
               if ~isempty(i3)
                  i2 = i3;
               end;
               i3=strfind(gmodel{i}(1:i2), '[');
               if isempty(i3)
                  gstat{ip1} = strtrim(gmodel{i}(1:i2 - 1));
               else
                  i5=i3+1;
                  i2=strfind(gmodel{i}(i5:i2), ']')+i5-1;
                  i4=strfind(gmodel{i}(i5:i2), ',')+i5-1;
                  while ~isempty(i4)
                     gstat{ip1}=strtrim(gmodel{i}(i5:i4 - 1));
                     ip1=ip1+1;
                     i5=i4+1;
                     i4=strfind(gmodel{i}(i5:i2), ',');
                  end;
                  gstat{ip1}=strtrim(gmodel{i}(i5:i2 - 1));
               end;  
               ip1 = ip1 + 1;
            end;
            end;
         end;
      end;
   end;
   [xx,ipars] = sort(pars);
   pars=pars(ipars);
   g_grind.solver.nonautonomous = 0;
%   g_grind.solver.isstochastic = 0;
%   for i = 1:size(gmodel, 2)
%      if ~isempty(strtrim(gmodel{i})) & (gmodel{i}(1) ~= '%')
%         i1 = strfind(gmodel{i}', rednoise(t,');
%         if isempty(i1)
%            i1 = strfind(gmodel{i}, 'randn(');
%            if isempty(i1)
%               i1 = strfind(gmodel{i}, 'rand(');
%           end;
%         end;
%         if ~isempty(i1)
%            if (i1 == 1) | isempty(strfind('1234567890abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_', gmodel{i}(i1 - 1)))
%               g_grind.solver.isstochastic = 1;
%            end;
%         end;
%      end;
%   end;
   i = 1;
   ip = 1;
   ip1 = ip1 - 1;
   while i <= size(pars, 1)
      p = char(pars(i));
      f = 0;
      if strcmp(p, 't')
         f = 1;
         g_grind.solver.nonautonomous = 1;
      end;
      for k = 1:ip1
         if strcmp(p, gstat{k})
            f = 1;
         end
      end
      if ~f
         g_grind.pars{ip} = p;
         ip = ip + 1;
      end
      i = i + 1;
      while i <= size(pars, 1) && strcmp(pars(i), p)
         i = i + 1;
      end
   end
   [xx,ipars] = sort(lower(g_grind.pars));
   g_grind.pars=g_grind.pars(ipars);
end;
