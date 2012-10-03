function s = i_disptext(s)
specials={'Delta','Gamma','Lambda','Omega','Phi','Pi','Psi','Sigma',...
   'Theta','Upsilon','Xi','alpha','beta','chi',...
   'delta','epsilon','eta','gamma','iota','kappa','lambda','mu','ni','nu',...
   'omega',  'phi','pi','psi','rho','sigma','tau','theta','upsilon',...
   'varsigma','vartheta','xi','zeta'};

nspec = 37; %length(specials);
syms=i_symvar(s,{'x','y','z','i','j'});

% suppress a warning if the string ends with '_'
f = strfind(s, '_');
if ~isempty(f) && (f(end)==length(s))
   s = [s(1:end - 1) '\_'];
end;

for j = 1:length(syms)
   i = 1;
   char = syms{j}(1);
   while (i < nspec) && (specials{i}(1) < char)
      i = i + 1;
   end;
   while (i<=nspec) && (specials{i}(1)==char)
      %   for i=1:length(specials)
      if strcmp(specials{i}, syms{j})
         s=strrep(s,specials{i},['{\' specials{i} '}']);
      end;
      i = i + 1;
   end;
   %  end;
end;
f = strfind(s, 'timesens(');
if ~isempty(f)
   f2=strfind(s, '''');
   s = [s(f2(1) + 1:f2(2) - 1) s(f2(2) + 2:end)];
end;
f = strfind(s, '_statevar(');
if ~isempty(f)
   s = [s(1:f(1)) ' of ' s(f(1) + 10:length(s) - 1)];
   s=strrep(s,'''','');
end;
f = strfind(s, 'outf(');
if ~isempty(f)
   f2 = strfind(s,',');
   if ~isempty(f2)
      fun = s(f(1) + 6:end);
      if strncmpi(fun, 'mov_', 4)
         fun = sprintf('moving %s', fun(5:end));
      elseif strncmpi(fun, 'perc', 4)
         fun = sprintf('%s%% percentile', fun(5:end));
      elseif strncmpi(fun, 'cover', 5)
         if (length(fun)>5)&&((fun(6)=='>')||(fun(6)=='<'))
            op = fun(6);
         else
            op = '>';
         end;
         if length(f2) == 2
            val = s(f2(2) + 1:end - 1);
            s = s(1:f2(2));
         else
            val = '0.01';
         end;
         fun = sprintf('cover%s%s', op, val);
      end;
      s = [fun ' of ' s(f2(1) + 1:length(s) - 1)];
   elseif strcmpi(s,'outf(''domeigen'')')
      s='dominant eigenvalue';
   elseif strcmpi(s,'outf(''realeigen'')')
      s='real eigenvalues';
   elseif strcmpi(s,'outf(''imageigen'')')
      s='imaginary eigenvalues';
   end;
   s=strrep(s,'''','');
end;

%this function always checks whether s2 is in s1


