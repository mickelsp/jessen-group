function i_sensplot(tag, data,parsets,vlabels,plabels)
if nargin == 0
   tag = 'init';
end;
switch tag
 case 'init'
   ud.output = data;
   ud.parsets = parsets;
   if exist('i_figno','file')
      n = i_figno('mcarlo') + 1;
   else
      n = 2;
   end;
   if ishandle(n)
      close(n);
   end;
   h0 = figure(n);
   set(h0,'Color',[0.914 0.914 0.914], ...
      'PaperPosition',[18 180 576 432], ...
      'PaperUnits','points', ...
      'Position',[269 251 560 420], ...
      'ResizeFcn','doresize(gcbf)', ...
      'Tag','Fig448', ...
      'ToolBar','figure', ...
      'DefaultaxesCreateFcn','plotedit(gcbf,''promoteoverlay''); ');
   set(h0, 'userdata', ud);
   uimenu('Parent',h0, ...
      'HandleVisibility','off', ...
      'Tag','ScribeHGBinObject', ...
      'Visible','off');
   uimenu('Parent',h0, ...
      'HandleVisibility','off', ...
      'Tag','ScribeFigObjStorage', ...
      'Visible','off');
   uimenu('Parent',h0, ...
      'HandleVisibility','off', ...
      'Tag','ScribeHGBinObject', ...
      'Visible','off');
   h1 = axes('Parent',h0, ...
      'Box','on', ...
      'CameraUpVector',[0 1 0], ...
      'Color',[1 1 1], ...
      'CreateFcn','', ...
      'Position',[0.1536 0.1762 0.7732 0.7024], ...
      'Tag','Axes1', ...
      'XColor',[0 0 0], ...
      'YColor',[0 0 0], ...
      'ZColor', [0 0 0]);
   line('Parent',h1, ...
      'Color',[0 0 1], ...
      'LineStyle','none', ...
      'Marker','.', ...
      'Tag','Axes1Line1', ...
      'XData',[], ...
      'YData', []);
   line('Parent',h1, ...
      'Color',[1 0 0], ...
      'LineStyle','-', ...
      'Marker','none', ...
      'Tag','line', ...
      'XData',[], ...
      'YData', []);
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'ListboxTop',0, ...
      'Callback','i_sensplot(''update'')', ...
      'Position',[101.25 9 99.75 15], ...
      'String',vlabels, ...
      'Style','popupmenu', ...
      'Tag','VarPopup', ...
      'Value', 1);
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[0.914 0.914 0.914], ...
      'ListboxTop',0, ...
      'Callback','i_sensplot(''update'')', ...
      'Position',[295.5 9 100.5 15], ...
      'String',plabels, ...
      'Style','popupmenu', ...
      'Tag','ParPopup', ...
      'Value', 1);
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[0.914 0.914 0.914], ...
      'ListboxTop',0, ...
      'Position',[239.25 6.75 45 15], ...
      'String','X-axis', ...
      'Style','text', ...
      'Tag','StaticText1');
   uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',[0.914 0.914 0.914], ...
      'ListboxTop',0, ...
      'Position',[47.25 8.25 45 15], ...
      'String','Y-axis', ...
      'Style','text', ...
      'Tag','StaticText1');
   i_sensplot('update')
 case 'update'
   h1=gcbf;
   if isempty(h1)
       h1=gcf;
   end;
   ud = get(h1, 'userdata');
   h=findobj(h1,'Tag','ParPopup');
   ipar = get(h, 'Value');
   pars = get(h, 'string');
   xlabel(pars{ipar});
   h=findobj(h1,'Tag','VarPopup');
   ivar = get(h, 'Value');
   vars = get(h, 'string');
   ylabel(vars{ivar});
   h=findobj(h1,'Tag','Axes1Line1');
   yy =  ud.output(:, ivar);
   xx = ud.parsets(:, ipar);
   set(h, 'XData', xx);
   set(h, 'YData', yy);
   xrange = [min(xx(~isnan(yy))) max(xx(~isnan(yy)))];
   yrange = [min(yy(~isnan(yy))) max(yy(~isnan(yy)))];
   if yrange(2) - yrange(1) > 0
      xlim(xrange);
      ylim(yrange);
      [sens,a,b]=i_mcsensanal(yy,xx);
      yregr = xrange' .* a + b;
      h=findobj(h1,'Tag','line');
      set(h, 'XData', xrange);
      set(h, 'YData', yregr);
   else
      sens = NaN;
   end;
  % h=findobj(h1,'Tag','TitleText');
   title(sprintf('Sensitivity "%s" for "%s" is %g',vars{ivar},strtrim(pars{ipar}),sens));
end

