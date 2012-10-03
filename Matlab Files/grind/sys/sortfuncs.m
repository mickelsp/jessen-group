function res=sortfuncs(funcs)
global g_grind;
if nargin==0
   funcs = str2cell(g_grind.funcs);
end;
if g_grind.statevars.vector
   filter = [g_grind.pars, g_grind.statevars.vectnames];
else
   filter = [g_grind.pars, g_grind.statevars.names];
end;
funs=cell(1,length(funcs));
for i = 1:length(funcs)
   f=strfind(funcs{i}, '=');
   if ~isempty(f)
      fun.comm = funcs{i};
      s = funcs{i};
      fun.name = strtrim(s(1:f(1) - 1));
      s = s(f(1) + 1:end);
      fun.references = i_symvar(s, filter);
      funs{i} = fun;
   end;
end;
%Sort
sortedfuncs=cell(1,length(funs));
refs=cell(1,length(funs));
j=0;
changedsome=1; %safety for errorous funcs;
while (j<length(funs)) && changedsome
   changedsome=0;
   for i = 1:length(funs)
      if ~isempty(funs{i}.name) && allunion(funs{i}.references,refs)
         changedsome=1;
         j=j+1;
         sortedfuncs{j}=funs{i}.comm;
         refs{j}=funs{i}.name;
         funs{i}.name='';
      end;
   end;
end;
res=sprintf('%s\n',sortedfuncs{:});
function res=allunion(s1,s2)
res=1;
for i=1:length(s1)
   found=0;
   j=1;
   while ~found && (j<=length(s2))
      if strcmp(s1{i},s2{j})
         found=1;
      end;
      j=j+1;
   end;
   if ~found
      res=0;
      return;
   end;
end;
         
      
      
