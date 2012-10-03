function [l_model, l_commands, l_scheme] = i_loadinifile(l_inifile)
global g_grind;
if ~isempty(strfind(l_inifile, '.all'))
   l_odefile=g_grind.odefile;
   load(l_inifile,'g_grind','-mat');
   l_model = g_grind.model;
   g_grind.odefile=l_odefile;
   l_commands = g_grind.commands;
else
   if isempty(strfind(l_inifile, '.'))
      l_inifile = [l_inifile '.ini'];
   end;
   fid = fopen(l_inifile, 'r');
   if fid > 0
      k = 1;
      model = 0;
      k2 = 1;
      k3=1;
      l_model = cell(1, 100);
      l_commands = cell(1, 100);
      l_scheme = cell(1, 100);
      while 1
         line = fgetl(fid);
         if ~ischar(line), break, end
         switch line
          case '%model'
            model = 1;
          case '%commands'
             model = 0;
         case '%scheme'
            model = 2;
          otherwise
            if ~isempty(line)
               if model==1
                  l_model{k} = line;
                  k = k + 1;
               elseif model==0
                  l_commands{k2} = line;
                  k2 = k2 + 1;
               else
                  l_scheme{k3} = line;
                  k3=k3+1;
               end
            end
         end;
      end;
      l_model = l_model(1:k - 1);
      l_commands = l_commands(1:k2 - 1);
      l_scheme = l_scheme(1:k3 - 1);
      fclose(fid);
   else
      l_model =[];
      l_commands =[];
      l_scheme=[];
   end;
end;
if nargout==0
   g_grind.model=l_model;
   g_grind.commands=l_commands;
   g_grind.inifile=l_inifile;
end;
