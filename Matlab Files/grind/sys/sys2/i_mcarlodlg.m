function  i_mcarlodlg(flag, par1)
evalin('base','global g_mcarlo;');
global g_grind g_mcarlo;
if nargin < 1
   flag = 'init';
end;
distrs={'Uniform','Normal','TruncNormal','LogNormal','paranal','hysteresis'};
switch flag
 case 'init'
   if ~isempty(g_mcarlo) && (~isfield(g_mcarlo, 'inifile') || strcmp(g_grind.inifile, g_mcarlo.inifile))
      ud = g_mcarlo.allpars;
      for i = 1:length(ud)
         ud{i}.value = evalin('base', ud{i}.name);
         if ~isnan(ud{i}.range)
            ud{i}.min = ud{i}.value .* (1 - ud{i}.range);
            ud{i}.max = ud{i}.value .* (1 + ud{i}.range);
         end;
      end;
   else
      g_mcarlo = [];
      ud = cell(size(g_grind.pars));
      for i = 1:length(g_grind.pars)
         p.name = g_grind.pars{i};
         p.descr = par('-d', p.name);
         if isempty(p.descr)
            p.descr = 'Parameter';
         end;
         p.value = evalin('base', p.name);
         p.selected = zeros(size(p.value));
         p.range = 0.1 + p.selected;
         p.min = p.value .* (1 - p.range);
         p.max = p.value .* (1 + p.range);
         p.sd = p.value .* p.range;
         p.distr = 'Uniform';
         ud{i} = p;
      end;
      npar = length(ud);
      if g_grind.statevars.vector
         for i = 1:length(g_grind.statevars.vectnames)
            p.name = g_grind.statevars.vectnames{i};
            p.descr = sprintf('Initial condition of %s', p.name);
            p.value = evalin('base', p.name);
            p.selected = zeros(size(p.value));
            p.range = 0.1 + p.selected;
            p.min = p.value .* (1 - p.range);
            p.max = p.value .* (1 + p.range);
            p.sd = p.value .* p.range;
            p.distr = 'Uniform';
            ud{npar + i} = p;
         end
      else
         for i = 1:g_grind.statevars.dim
            p.name = g_grind.statevars.names{i};
            p.descr = sprintf('Initial condition of %s', p.name);
            p.selected = 0;
            p.range = 0.1;
            p.value = evalin('base', p.name);
            p.min = p.value .* (1 - p.range);
            p.max = p.value .* (1 + p.range);
            p.sd = p.value .* p.range;
            p.distr = 'Uniform';
            ud{npar + i} = p;
         end;
      end
   end;
   pars = cell(1, length(ud));
   for i = 1:length(ud)
      pars{i} = ud{i}.name;
   end;
   
   h0 = figure('Color',[0.941 0.941 0.941], ...
      'MenuBar','none', ...
      'Name','Monte-Carlo sensitivity analysis', ...
      'NumberTitle','off', ...
      'PaperPosition',[18 180 576 432], ...
      'PaperUnits','points', ...
      'Position',[336 313 567 403], ...
      'Tag','Fig1', ...
      'ToolBar','none',...
      'CreateFcn',@(h,evnt)movegui(gcbf, 'center'));

    uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[1 1 1], ...
      'Callback','i_mcarlodlg(''listboxclick'')', ...
      'Position',[16.5 20.25 126 244.5], ...
      'String',pars, ...
      'Style','listbox', ...
      'Tag','ParList', ...
      'Value', 1);
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[0.941 0.941 0.941], ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[15 273.75 153 15], ...
      'String','Select parameters for sensitivity analysis:', ...
      'Style','text', ...
      'Tag','StaticText1');
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[0.941 0.941 0.941], ...
      'Callback','i_mcarlodlg(''selectedcheck'')', ...
      'ListboxTop',0, ...
      'Position',[157.5 225.75 72 21], ...
      'String','Selected?', ...
      'Style','checkbox', ...
      'Tag','SelectedCheck');
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[1 1 1], ...
      'Callback','i_mcarlodlg(''minedited'')', ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[258.75 142.5 93.75 18], ...
      'String','', ...
      'Style','edit', ...
      'Tag','MinEdit');
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[0.941 0.941 0.941], ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[157.5 141.75 65 15], ...
      'String','Minimum:', ...
      'Style','text', ...
      'Tag','MinText');
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[0.941 0.941 0.941], ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[157.5 105 65 15], ...
      'String','Maximum', ...
      'Style','text', ...
      'Tag','MaxText');
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[1 1 1], ...
      'Callback','i_mcarlodlg(''maxedited'')', ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[258.75 109.5 93.75 18], ...
      'String','', ...
      'Style','edit', ...
      'Tag','MaxEdit');
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[1 1 1], ...
      'Callback','i_mcarlodlg(''rangeedited'')', ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[258.75 78.75 93.75 18], ...
      'String','', ...
      'Style','edit', ...
      'Tag','RangeEdit');
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[0.941 0.941 0.941], ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[157.5 72 80 15], ...
      'String','Relative range:', ...
      'Style','text', ...
      'Tag','RangeText');
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[0.941 0.941 0.941], ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[158.25 260.25 200.25 15], ...
      'String','Parameter', ...
      'Style','text', ...
      'Tag','ParText');
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'Callback','i_mcarlodlg(''OKclick'')', ...
      'ListboxTop',0, ...
      'Position',[176.25 17.25 71.25 26.25], ...
      'String','OK', ...
      'Tag','OKButton');
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'Callback','i_mcarlodlg(''Cancelclick'')', ...
      'ListboxTop',0, ...
      'Position',[252 17.25 71.25 26.25], ...
      'String','Cancel', ...
      'Tag','CancelButton');
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'Callback','i_mcarlodlg(''Helpclick'')', ...
      'ListboxTop',0, ...
      'Position',[327.75 17.25 71.25 26.25], ...
      'String','Help', ...
      'Tag','HelpButton');
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'Callback','i_mcarlodlg(''SelectAllclick'')', ...
      'ListboxTop',0, ...
      'Position',[342.5 230.25 71.25 26.25], ...
      'String','Select all', ...
      'Tag','SelectallButtom');
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[1 1 1], ...
      'Callback','i_mcarlodlg(''distrclick'')', ...
      'ListboxTop',0, ...
      'Position',[258.75 175.5 91.5 18], ...
      'String',distrs, ...
      'Style','popupmenu', ...
      'Tag','DistrPopup', ...
      'Value', 1);
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[0.941 0.941 0.941], ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[157.5 175.5 45 12], ...
      'String','Distribution', ...
      'Style','text', ...
      'Tag','StaticText2');
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[1 1 1], ...
      'ListboxTop',0, ...
      'Position',[258.75 210.75 89.25 18], ...
      'String',' ', ...
      'Callback','i_mcarlodlg(''elemclick'')', ...
      'Style','popupmenu', ...
      'Visible','off',...
      'Tag','ElemList', ...
      'Value', 1);
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[0.941 0.941 0.941], ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[157.5 195.75 91.5 18], ...
      'String','Element vector/matrix', ...
      'Style','text', ...
      'Visible','off',...
      'Tag','ElemText');
   set(h0, 'userdata', ud);
   i_mcarlodlg('listboxclick', h0);
 case 'distrclick'
   ud = get(gcbf, 'userdata');
   h=findobj(gcbf,'tag','ParList');
   parnr = get(h, 'value');
   h=findobj(gcbf,'tag','DistrPopup');
   v = get(h, 'value');
   ud{parnr}.distr = distrs{v};
   set(gcbf, 'userdata', ud);
   i_mcarlodlg('listboxclick', gcbf);
 case 'listboxclick'
   if nargin >= 2
      h0 = par1;
   else
      h0 = gcbf;
   end;
   ud = get(h0, 'userdata');
   h=findobj(h0,'tag','ParList');
   parnr = get(h, 'value');
   if isempty(parnr)
      parnr = 1;
   end;
   h3=findobj(h0,'tag','ParText');
   set(h3,'string',sprintf('%s:  %s = %g',ud{parnr}.descr,ud{parnr}.name,mean(ud{parnr}.value(:))));
   h=findobj(h0,'tag', 'RangeEdit');
   set(h, 'string', num2str(ud{parnr}.range(1)));
   h=findobj(h0,'tag', 'MaxEdit');
   set(h, 'string', num2str(ud{parnr}.max(1)));
   h=findobj(h0,'tag','MinEdit');
   set(h, 'string', num2str(ud{parnr}.min(1)));
   h=findobj(h0,'tag', 'SelectedCheck');
   set(h, 'value', ud{parnr}.selected(1));
   h1=findobj(h0,'Tag','ElemList');
   h2=findobj(h0,'Tag','ElemText');
   if numel(ud{parnr}.value) > 1
      set(h3,'string',sprintf('%s:  mean(%s) = %g',ud{parnr}.descr,ud{parnr}.name,mean(ud{parnr}.value(:))));
      set(h1, 'visible', 'on');
      set(h2, 'visible', 'on');
      elems = cell(1, numel(ud{parnr}.value) + 1);
      if size(ud{parnr}.value, 2) > 1
         elems{1} = sprintf('%s(1:%d,1:%d)',ud{parnr}.name,size(ud{parnr}.value));
      else
         elems{1} = sprintf('%s(1:%d)', ud{parnr}.name, length(ud{parnr}.value));
      end;
      k = 2;
      if size(ud{parnr}.value, 2) > 1
         for i = 1:size(ud{parnr}.value, 1)
            for j = 1:size(ud{parnr}.value, 2)
               elems{k} = sprintf('%s(%d,%d)',ud{parnr}.name,i,j);
               k = k + 1;
            end;
         end;
      else
         for k = 1:size(ud{parnr}.value, 1)
            elems{k + 1} = sprintf('%s(%d)', ud{parnr}.name, k);
         end;
      end;
      set(h1, 'value', 1);
      set(h1, 'string', elems);
   else
      set(h1, 'visible', 'off');
      set(h1, 'value', 1);
      set(h2, 'visible', 'off');
   end;
   
   h=findobj(h0,'tag','DistrPopup');
   for i = 1:length(distrs)
      if strcmpi(distrs{i}, ud{parnr}.distr)
         set(h, 'value', i);
      end;
   end;
   if ~any(strcmpi({'uniform','paranal','hysteresis'},ud{parnr}.distr))
      h=findobj(h0,'tag', 'MaxEdit');
      set(h, 'visible', 'off');
      h=findobj(h0,'tag', 'MaxText');
      set(h, 'visible', 'off');
      h=findobj(h0,'tag', 'MinEdit');
      set(h, 'string', num2str(ud{parnr}.sd(1)));
      h=findobj(h0,'tag', 'MinText');
      set(h,'string','Stand.dev.');
      h=findobj(h0,'tag', 'RangeText');
      set(h,'string','CV (sd/mean)');
   elseif strcmpi('uniform', ud{parnr}.distr)
      h=findobj(h0,'tag', 'MaxEdit');
      set(h, 'visible', 'on');
      h=findobj(h0,'tag', 'MaxText');
      set(h,'string','Maximum:');
      set(h, 'visible', 'on');
      h=findobj(h0,'tag', 'MinText');
      set(h,'string','Minimum:');
      h=findobj(h0,'tag', 'RangeText');
      set(h,'string','Relative range:')
   else
      h=findobj(h0,'tag', 'MinText');
      set(h,'string','[Start, End]:');
       h=findobj(h0,'tag', 'MinEdit');
       if isfield(ud{parnr},'minmax')
          set(h,'string',sprintf('[%g %g]',ud{parnr}.minmax));
      else
          set(h,'string',sprintf('[%g %g]',ud{parnr}.min(1),ud{parnr}.max(1)));
          ud{parnr}.minmax=[ud{parnr}.min(1),ud{parnr}.max(1)];
      end;
      h=findobj(h0,'tag', 'MaxEdit');
      set(h, 'visible', 'on');
      if isfield(ud{parnr},'steps')
          set(h,'string',sprintf('%g',ud{parnr}.steps));
      else
          set(h,'string','40');
          ud{parnr}.steps=40;
      end;
      h=findobj(h0,'tag', 'MaxText');
      set(h,'string','Number of steps:');
      set(h, 'visible', 'on');
      h=findobj(h0,'tag', 'RangeText');
      set(h,'string','[Stabilizing, Writing]:')
      h=findobj(h0,'tag', 'RangeEdit');
      if isfield(ud{parnr},'nstabilwrite')
          set(h,'string',sprintf('[%g %g]',ud{parnr}.nstabilwrite));
      else
          set(h,'string','[600 300]');
          ud{parnr}.nstabilwrite=[600, 300];
      end;
      set(gcbf,'userdata',ud);
   end;
 case 'selectedcheck'
   ud = get(gcbf, 'userdata');
   h=findobj(gcbf,'tag','ParList');
   parnr = get(h, 'value');
   h=findobj(gcbf,'Tag','ElemList');
   elemnr =  get(h, 'value');
   h=findobj(gcbf,'tag', 'SelectedCheck');
   if elemnr == 1
      ud{parnr}.selected =  zeros(size(ud{parnr}.selected)) + get(h, 'value');
   else
      ud{parnr}.selected(elemnr - 1) =  get(h, 'value');
   end;
   set(gcbf, 'userdata', ud);
 case 'rangeedited'
   ud = get(gcbf, 'userdata');
   h=findobj(gcbf,'tag','ParList');
   parnr = get(h, 'value');
   h=findobj(gcbf,'Tag','ElemList');
   elemnr =  get(h, 'value');
   h=findobj(gcbf,'tag', 'RangeEdit');
   if any(strcmpi({'paranal','hysteresis'},ud{parnr}.distr))
      if elemnr  > 1
         error('GRIND:mcarlo:notimplemented','Vector parameters not yet implemented for paranal');
      else
         ud{parnr}.nstabilwrite = str2num(get(h, 'string')); %#ok
         set(gcbf, 'userdata', ud);
      end;
   else
      %if ~isnan(ud{parnr}.range)
      if elemnr == 1
         ud{parnr}.range = zeros(size(ud{parnr}.value)) + str2double(get(h, 'string'));
         ud{parnr}.min = ud{parnr}.value .* (1 - ud{parnr}.range);
         ud{parnr}.max = ud{parnr}.value .* (1 + ud{parnr}.range);
         ud{parnr}.sd = ud{parnr}.value .* ud{parnr}.range;
      else
         elemnr = elemnr - 1;
         ud{parnr}.range(elemnr) = str2double(get(h, 'string'));
         ud{parnr}.min(elemnr) = ud{parnr}.value(elemnr) .* (1 - ud{parnr}.range(elemnr));
         ud{parnr}.max(elemnr) = ud{parnr}.value(elemnr) .* (1 + ud{parnr}.range(elemnr));
         ud{parnr}.sd(elemnr) = ud{parnr}.value(elemnr) .* ud{parnr}.range(elemnr);
      end;
      %end;
      set(gcbf, 'userdata', ud);
      if strcmpi(ud{parnr}.distr, 'uniform')
         h=findobj(gcbf,'tag', 'MaxEdit');
         set(h, 'string', num2str(ud{parnr}.max(elemnr)));
         h=findobj(gcbf,'tag','MinEdit');
         set(h, 'string', num2str(ud{parnr}.min(elemnr)));
      else
         h=findobj(gcbf,'tag', 'MinEdit');
         set(h, 'string', num2str(ud{parnr}.sd(elemnr)));
      end;
   end
 case 'maxedited'
   ud = get(gcbf, 'userdata');
   h=findobj(gcbf,'tag','ParList');
   parnr = get(h, 'value');
   h=findobj(gcbf,'Tag','ElemList');
   elemnr =  get(h, 'value');
   h=findobj(gcbf,'tag', 'MaxEdit');
   if any(strcmpi({'paranal','hysteresis'},ud{parnr}.distr))
      if elemnr  > 1
         error('GRIND:mcarlo:notimplemented','Vector parameters not yet implemented for paranal');
      else
         ud{parnr}.steps = str2num(get(h, 'string')); %#ok
         set(gcbf, 'userdata', ud);
      end;
   else
      if elemnr == 1
         ud{parnr}.max = zeros(size(ud{parnr}.value)) + str2double(get(h, 'string'));
         if ud{parnr}.value == 0
            ud{parnr}.range = NaN;
         else
            r1 = iif(ud{parnr}.value == 0, 0, 1 - ud{parnr}.min ./ ud{parnr}.value);
            r2 = iif(ud{parnr}.value == 0, 0, ud{parnr}.max ./ ud{parnr}.value - 1);
            if r2 ~= r1
               ud{parnr}.range = NaN;
            else
               ud{parnr}.range = r1;
            end;
         end;
      else
         elemnr = elemnr - 1;
         ud{parnr}.max(elemnr) = str2double(get(h, 'string'));
         if ud{parnr}.value(elemnr) == 0
            ud{parnr}.range(elemnr) = NaN;
         else
            r1 = iif(ud{parnr}.value(elemnr) == 0, 0, 1 - ud{parnr}.min(elemnr) ./ ud{parnr}.value(elemnr));
            r2 = iif(ud{parnr}.value(elemnr) == 0, 0, ud{parnr}.max(elemnr) ./ ud{parnr}.value(elemnr) - 1);
            if r2 ~= r1
               ud{parnr}.range(elemnr) = NaN;
            else
               ud{parnr}.range(elemnr) = r1;
            end;
         end;
      end;
      set(gcbf, 'userdata', ud);
      h=findobj(gcbf,'tag', 'RangeEdit');
      set(h, 'string', num2str(ud{parnr}.range(elemnr)));
   end
 case 'minedited'
   ud = get(gcbf, 'userdata');
   h=findobj(gcbf,'tag','ParList');
   parnr = get(h, 'value');
   h=findobj(gcbf,'Tag','ElemList');
   elemnr =  get(h, 'value');
   h=findobj(gcbf,'tag', 'MinEdit');
   if any(strcmpi({'paranal','hysteresis'},ud{parnr}.distr))
      if elemnr  > 1
         error('GRIND:mcarlo:notimplemented','Vector parameters not yet implemented for paranal');
      else
         ud{parnr}.minmax = str2num(get(h, 'string')); %#ok
         set(gcbf, 'userdata', ud);
      end;
   else
      if elemnr == 1
         if ~strcmpi(ud{parnr}.distr, 'uniform')
            ud{parnr}.sd = zeros(size(ud{parnr}.value)) + str2double(get(h, 'string'));
            if ud{parnr}.value > 0
               ud{parnr}.range = iif(ud{parnr}.value == 0, 0, ud{parnr}.sd ./ ud{parnr}.value);
            else
               ud{parnr}.range = NaN;
            end;
         else
            ud{parnr}.min = zeros(size(ud{parnr}.value)) + str2double(get(h, 'string'));
            if ud{parnr}.value == 0
               ud{parnr}.range = NaN;
            else
               r1 = iif(ud{parnr}.value == 0, 0, 1 - ud{parnr}.min ./ ud{parnr}.value);
               r2 = iif(ud{parnr}.value == 0, 0, ud{parnr}.max ./ ud{parnr}.value - 1);
               if r2 ~= r1
                  ud{parnr}.range = NaN;
               else
                  ud{parnr}.range = r1;
               end;
            end;
         end;
      else
         elemnr = elemnr - 1;
         if ~strcmpi(ud{parnr}.distr, 'uniform')
            ud{parnr}.sd(elemnr) = str2double(get(h, 'string'));
            if ud{parnr}.value{elemnr} > 0
               ud{parnr}.range(elemnr) = iif(ud{parnr}.value(elemnr) == 0, 0, ud{parnr}.sd(elemnr) ./ ud{parnr}.value(elemnr));
            else
               ud{parnr}.range(elemnr) = NaN;
            end;
         else
            ud{parnr}.min(elemnr) = str2double(get(h, 'string'));
            if ud{parnr}.value{elemnr} > 0
               r1 = iif(ud{parnr}.value(elemnr) == 0, 0, 1 - ud{parnr}.min(elemnr) ./ ud{parnr}.value(elemnr));
               r2 = iif(ud{parnr}.value(elemnr) == 0, 0, ud{parnr}.max(elemnr) ./ ud{parnr}.value(elemnr) - 1);
               if r2 ~= r1
                  ud{parnr}.range(elemnr) = NaN;
               else
                  ud{parnr}.range(elemnr) = r1;
               end;
            else
               ud{parnr}.range(elemnr) = NaN;
            end;
         end;
      end;
      set(gcbf, 'userdata', ud);
      h=findobj(gcbf,'tag', 'RangeEdit');
      set(h, 'string', num2str(ud{parnr}.range(elemnr)));
   end;
 case 'elemclick'
   ud = get(gcbf, 'userdata');
   h1=findobj(gcbf,'tag','ParList');
   parnr = get(h1, 'value');
   h=findobj(gcbf,'tag','ElemList');
   elemnr = get(h, 'value');
   elems = get(h, 'string');
   h=findobj(gcbf,'tag','ParText');
   if elemnr == 1
      set(h,'string',sprintf('%s:  mean(%s) = %g',ud{parnr}.descr,ud{parnr}.name,mean(ud{parnr}.value(:))));
   else
      set(h,'string',sprintf('%s:  %s = %g',ud{parnr}.descr,elems{elemnr},ud{parnr}.value(elemnr-1)));
   end;
   if elemnr > 1
      elemnr = elemnr - 1;
   end;
   h=findobj(gcbf,'tag', 'RangeEdit');
   set(h, 'string', num2str(ud{parnr}.range(elemnr)));
   h=findobj(gcbf,'tag', 'MaxEdit');
   set(h, 'string', num2str(ud{parnr}.max(elemnr)));
   h=findobj(gcbf,'tag','MinEdit');
   if strcmpi(ud{parnr}.distr, 'uniform')
      set(h, 'string', num2str(ud{parnr}.min(elemnr)));
   else
      set(h, 'string', num2str(ud{parnr}.sd(elemnr)));
   end;
   h=findobj(gcbf,'tag', 'SelectedCheck');
   set(h, 'value', ud{parnr}.selected(elemnr));
 case 'SelectAllclick'
   ud = get(gcbf, 'userdata');
   Ans=inputdlg({'Rel.range/CV for all parameters (empty=no update):', ...
      'Distribution (Normal,Uniform,LogNormal,empty= no update)'},'Select all',1,{'0.1','Uniform'});
   if ~isempty(Ans)
      ran = str2double(Ans{1});
      if ~isempty(ran) && (ran > 0) && (ran < 1)
         for i = 1:length(ud)
            ud{i}.range = ran + zeros(size(ud{i}.value));
            ud{i}.min = ud{i}.value .* (1 - ud{i}.range);
            ud{i}.max = ud{i}.value .* (1 + ud{i}.range);
         end;
      end;
      distr = Ans{2};
      if ~isempty(distr)
         for i = 1:length(ud)
            ud{i}.distr = distr;
         end;
      end;
      for i = 1:length(ud)
         ud{i}.selected = ones(size(ud{i}.value));
      end;
      set(gcbf, 'userdata', ud);
      i_mcarlodlg('listboxclick');
   end;
 case 'OKclick'
   g_mcarlo.allpars = get(gcbf, 'userdata');
   g_mcarlo.pars = {};
   g_mcarlo.paranalpar =[];
   k = 1;
   nnegl = 0;
   nparanal =0;
   for j = 1:length(g_mcarlo.allpars)
      if g_mcarlo.allpars{j}.selected(1) && any(strcmpi({'paranal','hysteresis'},g_mcarlo.allpars{j}.distr))
         nparanal=nparanal+1;
         g_mcarlo.paranalpar=j;
      else
      for m = 1:size(g_mcarlo.allpars{j}.value, 1)
         for n = 1:size(g_mcarlo.allpars{j}.value, 2)
            if g_mcarlo.allpars{j}.selected(m, n)
               p = g_mcarlo.allpars{j};
               if  (size(g_mcarlo.allpars{j}.value, 2) > 1)
                  p.name = sprintf('%s(%d,%d)',p.name,m,n);
               elseif (numel(g_mcarlo.allpars{j}.value) > 1)
                  p.name = sprintf('%s(%d)', p.name, m);
               end;
               p.value = p.value(m, n);
               p.selected = p.selected(m, n);
               p.range = p.range(m, n);
               p.min = p.min(m, n);
               p.max = p.max(m, n);
               p.sd = p.sd(m, n);
               if ~isnan(p.sd) && (p.sd~=0)
                  g_mcarlo.pars{k} = p;
                  k = k + 1;
               else
                  nnegl = nnegl + 1;
               end;
            end;
         end;
      end;
      end;
   end;
   if nnegl > 0
      warning('GRIND:mcarlo:varsneglected','%d parameters/variables were neglected as their range was zero.\n', nnegl);
   end;
   if nparanal>1
       errordlg(sprintf('There are now %d parameters for paranal/hysteresis selected, select one maximal.',nparanal),'mcarlo')
   else
       close(gcbf);
   end;
 case 'Helpclick'
   disp('Help not yet implemented');
 case 'Cancelclick'
   close(gcbf);
end;
