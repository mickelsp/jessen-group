%SAVEPAR   Save current settings
%    Save the model, the current parameter settings and the initial
%    conditions.
%    Note that the current values of the parameters and the state variables
%    are the new defaults!
%   
%    Usage:
%    SAVEPAR - saves the model and parameter settings to the current ini file
%    SAVEPAR FILENAME - saves the model and settings to FILENAME
%    SAVEPAR FILENAME 1 - saves the model as FILENAME and overwrites
%    existing files without confirmation
%   
%    See also savemodel 

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function savepar(afile, overwrite)
global g_grind;
if ~isfield(g_grind, 'inifile')
   errordlg('No model selected');
   error('GRIND:savepar:NoModel','No model to save');
end;
if nargin == 0
   if isempty(g_grind.inifile)
      afile=uiputfile('*.ini;*.all', 'Save model and parameters as');
   else
      afile=uiputfile([g_grind.inifile ';*.all'], 'Save model and parameters as');
   end;
   overwrite = 1;
end;
if nargin  < 2
   overwrite = 0;
end;
overwrite = i_checkstr(overwrite);
if isempty(strfind(afile, '.all'))
   p = par('full');
   for j = 1:length(p)
      s=strtrim(p{j}(1:strfind(char(p{j}), '=') - 1));
      s2=strtrim(p{j}(strfind(char(p{j}), '=')  + 1:end - 1));
      if ~isempty(strfind(s2, '%'))
         s2=strtrim(s2(1:strfind(s2,'%')-1));
      end;
      f=strfind(s2, ';');
      if ~isempty(f)&&(f(1)==length(s2))
         s2=strtrim(s2(1:end-1));
      end;
      changecommand(s, p{j});
      changescheme(s,'%exp=', s2);
   end;
   [pg, comm]=par('groups');
   removecommands('par group ');
   if ~isempty(comm)
      g_grind.commands=[g_grind.commands , comm];
      for i=1:length(g_grind.pars)
          changescheme(g_grind.pars{i},'%tex=',g_grind.pargroups{i});
      end;
   end;
   N0 = i_initvar;
   if g_grind.statevars.vector
      for j = 1:length(g_grind.statevars.vectnames);
         kk = findcommand(g_grind.statevars.vectnames{j});
         if isempty(kk) || ~isempty(strfind(g_grind.commands{kk}, '['))
            ppar = evalin('base', g_grind.statevars.vectnames{j});
            s=sprintf('%s =[', g_grind.statevars.vectnames{j});
            for k = 1:size(ppar, 1)
               for l = 1:size(ppar, 2);
                  s = sprintf('%s%g, ',s,ppar(k,l));
               end;
               if j < size(ppar, 1)
                  s(length(s) - 1) = ';';
                  s = sprintf('%s...\n    ',s);
               end;
            end;
            s(length(s) - 1) = ']';
            s(length(s)) = ';';
            changecommand(g_grind.statevars.vectnames{j}, s);
         end;
      end
   else
      for j = 1:g_grind.statevars.dim
         s=sprintf('%s = %0.10g;', i_statevars_names(j), N0(j));
         changecommand(i_statevars_names(j), s);
         changescheme(i_statevars_names(j), '%exp=', sprintf('%0.10g', N0(j)));
      end;
   end;
end;
savemodel(afile, overwrite);

function icomm = findcommand(pname)
global g_grind;
icomm = [];
for k = 1:length(g_grind.commands)
   f=strfind(char(g_grind.commands{k}),'=');
   if ~isempty(f)
      s2 = strtrim(g_grind.commands{k}(1:f(1) - 1));
      if strcmp(pname, s2)
         icomm = k;
         return;
      end;
   end;
end;

function changecommand(pname, pnew)
global g_grind;
k = findcommand(pname);
if ~isempty(k)
   g_grind.commands{k} = pnew;
else
   %if not found then add
   g_grind.commands = [g_grind.commands, {pnew}];
end;

function removecommands(pname)
global g_grind;
g_grind.commands=g_grind.commands(~strncmp(pname,g_grind.commands,length(pname)));

function changescheme(pname, plabel, pnew)
global g_grind;
if isfield(g_grind, 'scheme')&&~strcmp(pnew,'[]')
   k = 1;
   while k < length(g_grind.scheme)
      s = g_grind.scheme{k};
      if strncmp(s, '%sym=', 5) && strcmp(s(6:end), pname)
         while (k<length(g_grind.scheme)) &&~(strncmp(g_grind.scheme{k},'%[',2) || strncmp(g_grind.scheme{k},plabel,5))
            k = k + 1;
         end;
         if (k < length(g_grind.scheme))  && strncmp(g_grind.scheme{k}, plabel, 5)
            s = g_grind.scheme{k};
            g_grind.scheme{k} = [s(1:5), pnew];
         end;
      end
      k = k + 1;
   end;
end;
