%VAL   initial and final values
%   Show initial and final values of state variables. Also to be used in a model
%   if you need to know the inital value of a state variable.
%   
%   Usage:
%   VAL  - displays the initial and final values of a state variable (after a run).
%   VAL('AVAR') - gets the inital value(s) of AVAR during a simulation.
%   VAL('AVAR(i,j)') - gets the inital value of AVAR(i,j) during a simulation.
%   
%   See also par, stat, model
%   

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function v=val(full, atitle)
global g_grind g_Y;
if nargin == 0
   full = '';
else
    if ~any(strcmp(full,{'full','?'}))
        %use in a model to get the initial value
        astatevar=full;
        v=evalin('base',astatevar);
        return;
    end;
end;
if ~isfield(g_grind,'statevars')||(g_grind.statevars.dim==0)
   disp('No state variables defined');
else
   disp(' ');
   n=g_grind.statevars.dim;
   if n>50&&~strcmp(full,'full')
      n=50;
      shortened=1;
   else
      shortened=0;
   end;
   if strcmp(full,'?')&&(n<20)
     vals=cell(n,1);
     N0 = i_initvar;
     for i=1:n
        vals{i}=num2str(N0(i));
     end;
     if nargin==1
        atitle='Set initial values';
     end;
     answer = inputdlg(g_grind.statevars.names, atitle, 1, vals);
     if ~isempty(answer)
        for i=1:n
           N0(i)=str2double(answer{i});
        end;
        i_keep(N0);
     end;
   end;    
   if isempty(g_Y)
      disp('Initial values of state variables:');
      N0 = i_initvar;
      maxparlen=par('-maxparlen');
      s=['%-' num2str(maxparlen) 's = %0.6g\n'];
       for i = 1:n
        fprintf(s,i_statevars_names(i),N0(i));
      end;
   else
      disp('Initial and final values of state variables in last run:');
      maxparlen=par('-maxparlen');
      s=['%-' num2str(maxparlen) 's = %0.5g / %0.5g\n'];
      for i = 1:n
         fprintf(s,i_statevars_names(i),g_Y(1, i),g_Y(end,i));
      end;
   end;
   if shortened
      fprintf('Only the first 50 of %d elements shown\n',g_grind.statevars.dim);      
   end;
end;
if isfield(g_grind,'permanent')&&~isempty(g_grind.permanent)
   disp(' ');
  defpermanent('-l');
end;
