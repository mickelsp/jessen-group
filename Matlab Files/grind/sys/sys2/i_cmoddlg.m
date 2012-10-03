% actions for the MODEL dialog
function res = i_cmoddlg(action)
global g_grind;
f = gcf;
switch action
 case 'makemodel'
   if ~isempty(gcbf)
      h = get(gcbf, 'Userdata');
      if ~isempty(h) && isfield(h, 'clear') && h.clear
         i_startup;
         initgrind;
         evalin('base','clear global;');
         h = get(gcbf, 'Userdata');
      else
         h.check = 1;
      end;
   else
      h.check  = 1;
   end;
   if h.debug
      l = i_cmoddlg('getmodel');
      i_mmodel(l);
      i_parcheck(1);
   else
      try
         err = 0;
         s = '';
         try
            l = i_cmoddlg('getmodel');
            i_mmodel(l);
         catch %#ok
            try
               [err, s] = i_parcheck(1);
               err = ~err;
            catch err1
        %       err1 = lasterror;
               s = striphtml(err1.message);
               err = 1;
            end;
         end;
         if ~err && h.check
            [err, s] = i_parcheck(1);
            err = ~err;
         end;
         if err && h.check
            if i_myerrordlg(s, 'Error in model') == 2
               close(gcf);
               i_startup;
               if ~isfield(h, 'era')
                  i_use(l.inifile, 1, 1, 0, h.era)
               else
                  i_use(l.inifile, 1, 1, 0, 1);
               end;
            else
               model(l.inifile);
               return;
            end;
         end;
      catch err
         %err = lasterror;
         %    errordlg(striphtml([s err.message]));
          rethrow(err);
      end;
   end;
 case 'undoclearScheme'
   if isfield(g_grind, 'oldscheme')
      g_grind.scheme = g_grind.oldscheme;
      g_grind = rmfield(g_grind, 'oldscheme');
      h = get (gcbf, 'Userdata');
      h.scheme = g_grind.scheme;
      set(gcbf, 'Userdata', h);
      h=findobj(gcbf,'Tag','ClearSchemeButton');
      set(h,'String','Clear scheme');
      set(h,'TooltipString','Clear the vismod scheme');
      set(h,'Callback','i_cmoddlg(''clearScheme'')');
      h=findobj(gcbf,'Tag','modellabel');
      set(h,'String','Model equations (use vismod to change)');
      h=findobj(gcbf,'Tag','Model');
      set(h,'Enable','off');
   end;
 case 'clearScheme'
   button = questdlg('Are you sure you want to delete the scheme?', ...
      'model','Yes', 'No','No');
   if strcmp(button, 'Yes')
      g_grind.oldscheme = g_grind.scheme;
      h = get (gcbf, 'Userdata');
      set(gcbf,'Userdata',rmfield(h,'scheme'));
      h=findobj(gcbf,'Tag','ClearSchemeButton');
      set(h,'String','Undo clear');
      set(h,'TooltipString','Undo clearing of the vismod scheme');
      set(h,'Callback','i_cmoddlg(''undoclearScheme'')');
      h=findobj(gcbf,'Tag','modellabel');
      set(h,'String','Model equations (diffential/difference equations) ');
      h=findobj(gcbf,'Tag','Model');
      set(h,'Enable','on');
   end;
 case 'resize'
   f = gcbf;
   figwidth = 430;
   figheight = 320;
   p = get(f, 'position');
   positionobj(f, 'OpenFileButton', [0 1 1 0], [339 276 83 25], [figwidth figheight], p(3:4));
   positionobj(f, 'namelabel', [1 0 1 0], [11 281 106 16], [figwidth figheight], p(3:4));
   positionobj(f, 'FileName', [1 1 1 0], [130 284 181 16], [figwidth figheight], p(3:4));
   positionobj(f, 'CDButton', [0 1 1 0], [314 284 14 16], [figwidth figheight], p(3:4));
   positionobj(f, 'modellabel', [1 0 1 0], [11 255 295 16], [figwidth figheight], p(3:4));
   positionobj(f,'ClearSchemeButton', [1 0 1 0], [200 260 70 16], [figwidth figheight], p(3:4));
   posmod = positionobj(f, 'Model', [1 1 1 1], [11 143 317 116], [figwidth figheight], p(3:4));
   posparlab = positionobj(f, 'parameterslabel', [1 0 0 1], [11 120 295 16], [figwidth figheight], p(3:4));
   pospar = positionobj(f, 'Parameters', [1 1 0 1], [11 6 318 118], [figwidth figheight], p(3:4));
   positionobj(f, 'ExtractButton', [0 1 0 1], [339 145 83 25], [figwidth figheight], p(3:4));
   positionobj(f, 'CancelButton', [0 1 0 1], [339 46 83 25], [figwidth figheight], p(3:4));
   positionobj(f, 'OKButton', [0 1 0 1], [339 13 83 25], [figwidth figheight], p(3:4));
   positionobj(f, 'OpenFileButton', [0 1 1 0], [339 276 83 25], [figwidth figheight], p(3:4));
   positionobj(f, 'NewButton', [0 1 0 1], [339 111  83 25], [figwidth figheight], p(3:4));
   positionobj(f, 'Helpbutton', [0 1 0 1], [339 78 83 25], [figwidth figheight], p(3:4));
   if ~isempty(posmod) && (posmod(4) ~= pospar(4))
      height = (posmod(4) + pospar(4)) / 2;
      %    if height>posmod(4)
      %      trans=posmod(2)-(posmod(4)-height);
      %      posmod(2)=posmod(2)+trans;
      %   else
      trans = height - pospar(4);
      posmod(2) = posmod(2) + trans;
      posparlab(2) = posparlab(2) + trans;
      %pospar(2)=pospar(2)-trans;
      %  end;
      pospar(4) = height;
      posmod(4) = height;
      set(findobj(f,'tag','Model'),'position',posmod);
      set(findobj(f,'tag','Parameters'),'position',pospar);
      set(findobj(f,'tag','parameterslabel'),'position',posparlab);
   end;
 case 'init'
   f = gcf;
   ud.clear = 0;
   ud.updated = 0;
   ud.era = 1;
   ud.scheme = {};
   set(f, 'Userdata', ud);
   if ~isempty(g_grind) && isfield(g_grind, 'model') && ~isempty(g_grind.model)
      %      global g_grind;
      if isfield(g_grind, 'scheme')  && ~isempty(g_grind.scheme)
         ud.scheme = g_grind.scheme;
         set(f, 'Userdata', ud);
      end;
      h = findobj(f,'Tag','FileName');
      set(h, 'String', g_grind.inifile);
      h = findobj(f,'Tag','Model');
      for i = 1:length(g_grind.model)
         if ~isempty(g_grind.model{i})
            g_grind.model{i} = strtrim(g_grind.model{i});
         end;
         if isempty(g_grind.model{i})
            g_grind.model{i} = ' ';
         else
         end;
      end;
      set(h, 'String', g_grind.model);
      h = findobj(f,'Tag','Parameters');
      for i = 1:length(g_grind.commands)
         if isempty(g_grind.commands{i})
            g_grind.commands{i} = ' ';
         end;
      end;
      set(h, 'String', g_grind.commands);
   end;
 case 'getmodel'
   f = gcbf;
   h = findobj(f, 'Tag', 'Model');
   r.model =  i_memo2cell(get(h, 'String'));
   h = findobj(f,'Tag','Parameters');
   r.commands = i_memo2cell(get(h, 'String'));
   h = findobj(f, 'Tag', 'FileName');
   r.inifile = strtrim(get(h, 'String'));
   h = get(f, 'Userdata');
   r.scheme = {};
   if ~isempty(h) && isfield(h, 'scheme')
      r.scheme = h.scheme;
   end;
   if ~isempty(h) && isfield(h, 'updated') && h.updated
      if isempty(r.inifile)
         s = inputdlg({'Enter file name'}, 'File name', 1, {'curr_mod.ini'});
         r.inifile = s{1};
      end;
      g_grind.model = r.model;
      g_grind.commands = r.commands;
      g_grind.inifile = r.inifile;
      g_grind.scheme = r.scheme;
      savemodel(char(r.inifile));
   end;
   if nargout == 1
      res = r;
   else
      g_grind.model = r.model;
      g_grind.commands = r.commands;
      g_grind.inifile = r.inifile;
      g_grind.scheme = r.scheme;
   end;
 case 'newmodel'
   savechanges;
   setupdated(0);
   s = ' ';
   h = findobj(f, 'Tag', 'Model');
   set(h, 'String', s);
   h = findobj(f, 'Tag', 'Parameters');
   set(h, 'String', s);
   h = findobj(f, 'Tag', 'FileName');
   set(h, 'String', s);
 case 'updated'
   setupdated(1)
 case 'extract'
   %   global g_grind;
   if ~isempty(g_grind) && isfield(g_grind, 'model')
      oldgmodel = g_grind.model;
      if isfield(g_grind, 'pars')
         oldgpars = g_grind.pars;
      else
         oldgpars = [];
      end;
   else
      oldgmodel = [];
      oldgpars = [];
   end;
   g_grind.pars = [];
   h = findobj(f, 'Tag', 'Model');
   g_grind.model = i_memo2cell(get(h, 'String'));
   i_fillpars(1);
   h = findobj(f, 'Tag', 'Parameters');
   s1 = i_memo2cell(get(h, 'String'));
   hh = length(s1);
   while (hh > 0) && isempty(strtrim(s1{hh}))
      hh = hh - 1;
   end;
   s = cell(1,size(g_grind.pars, 2) + 10+hh);
   s2 = cell(hh, 1);
   k1=0;
   for i = 1:hh
      f1 = strfind(s1{i},'=');
      if ~isempty(f1)
         k1 = k1 + 1;
         s2{k1} = strtrim(s1{i}(1:f1(1) - 1));
      end;
      s{i} = s1{i};
   end;
   s2 = s2(1:k1);
   done = 0;
   for i = 1:size(g_grind.pars, 2)
      f1 = 0;
      for j = 1:length(s2)
         if strcmp(g_grind.pars{i}, s2{j})
            f1 = 1;
            hh = hh - 1;
         end;
      end;
      if ~f1
         done = 1;
         s{i + hh} = [char(g_grind.pars{i}) '='];
      end;
   end;
   j = length(s);
   while (j > 0) &&  isempty(s{j})
      j = j - 1;
   end;
   s = s(1:j);
   set(h, 'String', s);
   g_grind.model = oldgmodel;
   g_grind.pars = oldgpars;
   if ~done
      msgbox('No new parameters found');
   else
      setupdated(1)
   end;
 case 'keydown'
   ch = get(f, 'CurrentCharacter');
   if ~isempty(ch)
      if (ch == char(4)) %control - d sets in debug mode
         ud = get(f, 'userdata');
         ud.debug = ~ud.debug;
         set(f, 'userdata', ud);
         if ud.debug
            msgbox('Ctrl-D pressed: Debug mode','model')
         else
            msgbox('Ctrl-D pressed: Normal mode: error messages displayed','model');
         end;
      elseif ch == char(8) %ctrl - h
         i_cmoddlg('help');
         
      elseif ch == char(3) %ctrl - c
         if strcmp('Yes', questdlg('OK to cancel?','Model','Yes','No','No'))
            i_cmoddlg('cancel');
         end;
      elseif ch == char(13) %enter
         i_cmoddlg('ok');
      end;
   end;
 case 'ok'
   set(f, 'DeleteFcn', 'i_cmoddlg(''makemodel'')' );
   delete(f);
 case 'cancel'
   if isempty(g_grind) || (isfield(g_grind, 'model') && isempty(g_grind.model))
      clear g_grind;
   end;
   delete(f);
 case 'cd'
   setupdated(1)
   h=findobj(f,'Tag','FileName');
   curfil = get(h, 'String');
   if isempty(curfil)
      curfil = '*.ini';
   end;
   if isempty(strfind(curfil, '.ini'))
      curfil = [curfil '.ini'];
   end;
   [filename, pathname] = uiputfile(curfil,'Select file name to save model');
   if ~isempty(filename) && ischar(filename)
      h=findobj(f,'Tag','FileName');
      set(h, 'String', char(filename));
   end;
   if ~isempty(pathname) && ischar(pathname)
      cd(pathname);
   end;
 case 'openfile'
   savechanges;
   [filename, pathname] = uigetfile('*.ini','Select file name');
   h=findobj(f,'Tag','FileName');
   if (filename ~= 0)
      set(h, 'String', char(filename));
      cd(pathname);
      i_cmoddlg('loadfile');
   end
 case 'loadfile'
   setupdated(0);
   h = findobj(f, 'Tag', 'FileName');
   l_inifile = strtrim(get(h, 'String'));
   [l_model, l_commands, l_scheme] = i_loadinifile(l_inifile);
   if ~isempty(l_model)
      h = findobj(f, 'Tag','Model');
      set(h, 'String', l_model);
      h = findobj(f, 'Tag', 'Parameters');
      set(h, 'String', l_commands);
      ud = get(f, 'Userdata');
      ud.scheme = l_scheme;
      set(f, 'Userdata', ud);
   else
      h = findobj(f, 'Tag', 'Model');
      set(h, 'String', {'', ''});
      h = findobj(f, 'Tag', 'Parameters');
      set(h, 'String', {'', ''});
      ud = get(f, 'Userdata');
      ud.scheme = {};
      set(f, 'Userdata', ud);
   end;
end;

function setupdated(value)
f = gcf;
ud = get(f, 'Userdata');
ud.updated = value;
set(f, 'Userdata', ud);
return;

function savechanges
global g_grind;
f = gcf;
ud = get(f, 'Userdata');
if ud.updated
   h = findobj(f, 'Tag', 'FileName');
   l_inifile = get(h, 'String');
   if ~isempty(l_inifile)&& strcmp('Yes',questdlg(['Do you want to save changes in ' l_inifile '?'], ...
         'Save changes','Yes','No','Yes'))
      oldgmodel = g_grind.model;
      oldgcommands = g_grind.commands;
      oldginifile = g_grind.inifile;
      h = findobj(f, 'Tag', 'Model');
      g_grind.model = transpose(get(h, 'String'));
      h=findobj(f, 'Tag', 'Parameters');
      g_grind.commands = transpose(get(h, 'String'));
      savemodel(l_inifile, 1);
      g_grind.model = oldgmodel;
      g_grind.commands = oldgcommands;
      g_grind.inifile = oldginifile;
   end;
end;
function pos = positionobj(f, tag, adj, defsize, oldsize, newsize)
% adj = left right top bottom
% position = left bottom width height
h = findobj(f, 'tag', tag);
if ~isempty(h)
   pos = defsize;
   if adj(2) && ~adj(1)
      pos(1) = newsize(1) - (oldsize(1) - pos(1));
   end;
   if adj(2) && adj(1)
      pos(3) = pos(3) - (oldsize(1) - newsize(1));
      if pos(3) < 1
         pos(3) = 1;
      end;
   end;
   if adj(3) && ~adj(4)
      pos(2) = newsize(2) - (oldsize(2) - pos(2));
   end;
   if adj(4) && adj(3)
      pos(4) = pos(4) - (oldsize(2) - newsize(2));
      if pos(4) < 1
         pos(4) = 1;
      end;
   end;
   set(h, 'position', pos);
else
   pos = [];
end;


