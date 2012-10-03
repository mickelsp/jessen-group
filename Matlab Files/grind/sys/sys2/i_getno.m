function iX = i_getno(avar)
global g_grind;
[iX.no, iX.vecno] = i_varno(avar);
iX.isvar = ~isempty(iX.no);
iX.ispar = 0;
iX.isfun = 0;
iX.isext=0;
iX.isperm = 0;
iX.ndx = [];
if ~iX.isvar
   iX.no = i_parno(avar);
   if ~isempty(iX.no)
      iX.ispar = 1;
   else
   %   if ~isfield(g_grind.funcnames, 'dims')
   %      i_evalfuncs;
   %   end;
      iX.no = i_permno(avar);
      iX.isperm = ~isempty(iX.no);
      if ~iX.isperm
         [iX.no, iX.vecno] = i_funno(avar);
         iX.isfun = ~isempty(iX.no);
         if ~iX.isfun
            iX.no = i_externno(avar);
            iX.isext = ~isempty(iX.no);
         end;
      end;
   end;
end;
if ~iX.isfun
   f =  strfind(avar, '(');
   s='[]';
   if ~isempty(f)
      f2 = strfind(avar, ')');
      f1 = strfind(avar, ',');
      if isempty(f1)
         iX.ndx = str2double(avar(f + 1:f2 - 1));
      else
         if iX.ispar || iX.isvar || iX.isperm
            s = sprintf('sub2ind(size(%s),%s',avar(1:f - 1),avar(f + 1:end));
         else
            if ~isfield(g_grind.funcnames, 'dims')
               s='NaN';
            elseif ~isempty(iX.no)
               s = sprintf('sub2ind([%d,%d],%s',g_grind.funcnames.dims{iX.no}.dim1,g_grind.funcnames.dims{iX.no}.dim2,avar(f + 1:end));
            end;
         end;
         iX.ndx = evalin('base', s);
      end;
   end;
   if iX.isfun
      iX.no = iX.no + iX.ndx - 1;
   end;
end;
function permno = i_permno(apar)
global g_grind;
permno = [];
f =  strfind(apar, '(');
if ~isempty(f)
   apar = apar(1:f(1) - 1);
end;
if isfield(g_grind,'permanent')
for k = 1:size(g_grind.permanent, 2)
   if strcmp(apar, char(g_grind.permanent{k}.name))
      permno = k;
      return;
   end
end
end;
function extno = i_externno(apar)
global g_grind;
extno = [];
f =  strfind(apar, '(');
if ~isempty(f)
   apar = apar(1:f(1) - 1);
end;
for k = 1:size(g_grind.externvars, 2)
   if strcmp(apar, char(g_grind.externvars{k}.name))
      extno = k;
      return;
   end
end

function parno = i_parno(apar)
global g_grind;
parno = [];
f =  strfind(apar, '(');
if ~isempty(f)
   apar = apar(1:f(1) - 1);
end;
for k = 1:length(g_grind.pars)
   if strcmp(apar, char(g_grind.pars{k}))
      parno = k;
      return;
   end
end
function [no, funno] = i_funno(apar)
global g_grind;
funno = [];
no = [];
f =  strfind(apar, '(');
f2 = strfind(apar, ')');
if (length(f)==1) && (length(f2)==1) && (f2==length(apar))
   apar = apar(1:f(1) - 1);
elseif ~isempty(f)
   return;
end;
for k = 1:size(g_grind.funcnames.names, 2)
   if strcmp(apar, char(g_grind.funcnames.names{k}))
      funno = k;
      if ~isfield(g_grind.funcnames, 'dims')
         no=NaN;
      else
         no = g_grind.funcnames.dims{k}.from;
      end;
      return;
   end
end
