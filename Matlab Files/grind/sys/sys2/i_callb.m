function i_callb(flag, avar, ndays)
global g_t g_Y g_grind;
switch flag
 case 'mdown'
   %  ax = get(gcbf, 'CurrentAxes');
   N0 = i_initvar;
   old_Y = g_Y;
   oldN0 = N0;
   ax = get(gcbf, 'CurrentAxes');
   initpt = get(ax, 'CurrentPoint');
   XLim = get(ax, 'Xlim');
   YLim = get(ax, 'Ylim');
   Ylabel= get(get(ax,'ylabel'),'string');
   Zlabel= get(get(ax,'zlabel'),'string');

   labels={get(get(gca,'xlabel'),'string'),get(get(gca,'ylabel'),'string'),get(get(gca,'zlabel'),'string')};
   if ~isempty(labels{3})
      initpt=i_depthslider('init',initpt,labels);
      if isempty(initpt)
         return;
      end;
   else
      initpt = initpt(1, [1, 2]);
   end;
   p2 = [];
   ndx = 1:g_grind.statevars.dim;
   if ~((initpt(1) < XLim(1)) || (initpt(1) > XLim(2)) || (initpt(2) < YLim(1)) || (initpt(2) > YLim(2)))
      %g_Y = transpose(N0);
      ix = i_getno(g_grind.xaxis.var);
      p1 = i_getparvar(ix);
      N0 = i_setparvar(ix, N0, initpt(1));
      ndx=ndx(ndx ~= ix.no);
      iy = i_getno(g_grind.yaxis.var);
      if isempty(strfind(char(Ylabel),''''))
         p2 = i_getparvar(iy);
         N0 = i_setparvar(iy, N0, initpt(2));
         ndx=ndx(ndx ~= iy.no);
      end
      if length(initpt)>2
            iz = i_getno(g_grind.zaxis.var);
        if isempty(strfind(char(Zlabel),''''))
           p2 = i_getparvar(iz);
           N0 = i_setparvar(iz, N0, initpt(3));
           ndx=ndx(ndx ~= iz.no);
        end
      end;
      mess = cell(1, g_grind.statevars.dim + 1);
      mess{1} = 'Initial conditions set at mouse cursor:';
      i_keep(transpose(N0))
      for i = 1:g_grind.statevars.dim
         mess{i + 1}=[i_statevars_names(i) ' = ' num2str(N0(i)) ];
      end;
      i = i + 1;
      if ix.ispar
         i = i + 1;
         mess{i}=[g_grind.pars{ix.no} ' = ' num2str(initpt(1))];
      end;
      if iy.ispar
         i = i + 1;
         mess{i}=[g_grind.pars{iy.no} ' = ' num2str(initpt(2))];
      end;
      H = gcbf;
      ud = get(H, 'userdata');
      if isfield(ud, 'stop')
         ud.stop = 0;
         set(H, 'UserData', ud);
         set(findobj(H,'Tag','stop'),'Visible','off');
         drawnow;
      end;
      set(H, 'WindowButtonDown', '');
      set(H, 'WindowButtonMotionFcn', '');
      oldpoint = get(H, 'Pointer');
      set(H,'Pointer','arrow');
      try
         opt = i_mymenu(mess);
      catch %#ok
         opt = 0;
         resetpoint(H, oldpoint);
      end;
      resetpoint(H, oldpoint);
      if isfield(ud, 'iters')
         g_grind.solver.iters = g_grind.solver.iters * ud.iters;
      end;
      % try
      if ~isempty(ndx) && ((opt==1) || (opt==2))
         vals =  cell(1, length(ndx));
         for i = 1:length(ndx)
            vals{i} = num2str(N0(ndx(i)));
         end;
         answer = inputdlg({i_statevars_names(ndx)}, 'Not all initial conditions are on the axes, if needed edit these', 1, vals);
         if ~isempty(answer)
            for i = 1:length(answer)
               N0(ndx(i)) = str2double(answer{i});
            end;
         end;
         i_keep(transpose(N0));
      end;
      switch opt
       case 1
         if isfield(g_grind, 'tilman')
            BtnName=questdlg('Do you want to set the resource supply parameters to the same value?','Yes');
            if strcmp(BtnName, 'Yes')
               try
                  evalin('base', g_grind.tilman.setS);
               catch err
  %                err=lasterror;
                  disp('error in evaluating g_grind.tilman.setS');
                  rethrow(err);
               end;
            end;
            ru;
         else
            ru;
         end;
       case 2
         findeq;
       case 4
         docancel(ix, iy, oldN0, p1, p2, old_Y);
         plotedit on;
         if ishandle(H)
            set(H,'Pointer','arrow');
         end;
       case 5
         docancel(ix, iy, oldN0, p1, p2, old_Y);
      end
      if isfield(ud, 'iters')
         g_grind.solver.iters = g_grind.solver.iters / ud.iters;
      end;
      %  catch
      %     if isfield(ud, 'iters')
      %        g_grind.solver.iters = g_grind.solver.iters / ud.iters;
      %    end;
      % end
   end;
 case 'mmove'
   f = gcbf;
   ax = get(f, 'CurrentAxes');
   initpt = get(ax, 'CurrentPoint');
   Lims = zeros(3, 2);
   Lims(1,:) = get(ax, 'Xlim');
   Lims(2,:) = get(ax, 'Ylim');
   Lims(3,:) = get(ax, 'Zlim');
   axlabels={'xlabel','ylabel','zlabel'};
   fax=find(abs(initpt(1, :) - initpt(2, :)) <= 0.05);
   nval = 0;
   if length(fax) > 1
      s = sprintf(' | ');
      for i = 1:length(fax)
         if ~((initpt(1, fax(i)) < Lims(fax(i), 1)) || (initpt(1, fax(i)) > Lims(fax(i), 2)))
            lab = get(get(ax, axlabels{fax(i)}), 'string');
            lab = clearcode(lab);
            if isempty(lab) || (length(lab) > 8)
               lab = [axlabels{fax(i)}(1) '-axis'];
            end;
            nval = nval + 1;
            s=sprintf('%s  %s = %0.3g', s, lab, initpt(1, fax(i)));
         end;
      end;
   else
      s = sprintf(' | ');
      fax=1:3;
      for i = 1:length(fax)
         if ~((initpt(1, fax(i)) < Lims(fax(i), 1)) || (initpt(1, fax(i)) > Lims(fax(i), 2)))
            lab = get(get(ax, axlabels{fax(i)}), 'string');
            lab = clearcode(lab);
            if isempty(lab) || (length(lab) > 8)
               lab = [axlabels{fax(i)}(1) '-axis'];
            end;
            nval = nval + 1;
            s=sprintf('%s  %s = %0.3g', s, lab, initpt(1, fax(i)));
         end;
      end;
   end;
   if nval > 1
      set(f,'Pointer','crosshair')
   else
      s = '';
      set(f,'Pointer','arrow');
   end;
   t = get(f, 'name');
   f1 = strfind(t, ' | ');
   if ~isempty(f1)
      t = t(1:f1(1) - 1);
   end;
   set(f, 'name', [t s]);
 case 'mdown2'
   %  global g_t g_Y g_grind;
   N0 = i_initvar;
   ax = get(gcbf, 'CurrentAxes');
   initpt = get(ax, 'CurrentPoint');
   initpt = initpt(1, [1, 2]);
   XLim = get(ax, 'Xlim');
   YLim = get(ax, 'Ylim');
   if ~((initpt(1) < XLim(1)) || (initpt(1) > XLim(2)) || (initpt(2) < YLim(1)) || (initpt(2) > YLim(2)))
      tclick = initpt(1);
      iy = i_getno(avar);
      if ~isempty(iy.no)
         N0(iy.no) = initpt(2);
      end;
      i_ru(g_grind.odefile, tclick, ndays - tclick, N0, 1);
      %[g_t,g_Y]=feval(g_grind.solver.name,g_grind.odefile,[tclick; t+ndays],N0,g_grind.solver.opt);
      hold on;
      plot(g_t, g_Y(:, iy.no), g_grind.pen.pen, 'Color', g_grind.pen.color);
      if g_grind.statevars.dim > 1
         title(['Other intial conditions ' i_othervars(N0, iy.no)]);
      end;
      nextpen;
   end;
 case 'keypressed'
   key = double(get(gcbf, 'CurrentCharacter'));
   if isempty(key)
      return
   end
   switch key
    case 1 %Ctrl - A
      ax = get(gcbf, 'CurrentAxes');
      initpt = get(ax, 'CurrentPoint');
      if ~isempty(initpt)
         i_draw2Darrow(initpt(1, 1), initpt(1, 2));
      end;
    case 65 %Shift - A
      checkdlg;
      ax = get(gcbf, 'CurrentAxes');
      initpt = get(ax, 'CurrentPoint');
      if ~isempty(initpt)
         arrows('delnearest', [initpt(1, 1), initpt(1, 2)]);
     end;
    case 2 %Ctrl - B
      checkdlg;
      backw;
    case 11 %Ctrl - K
      checkdlg;
      ke;
    case 18 %Ctrl - R
      checkdlg;
      ru;
    case 20 %Ctrl - T
      checkdlg;
      time;
    case 21 %Ctrl - U
      conteq('-r');
      figure(i_figno('phase2'));
    case 27 %Esc
      f = gcbf;
      if f == i_figno('dialog')
         close(f)
      else
         ud = get(f, 'userdata');
         if ~isempty(ud) && isfield(ud, 'stop');
            ud.stop = 1;
            set(f, 'UserData', ud);
         end;
      end;
    case 5 %Ctrl - E
      checkdlg;
      eigen;
    case 64 %Shift - 2
      checkdlg;
      null;
    case 35 %Shift - 3
      checkdlg;
      null3;
    case 25 %Ctrl - Y
      checkdlg;
      ax = get(gcbf, 'CurrentAxes');
      initpt = get(ax, 'CurrentPoint');
      if ~isempty(initpt)
         arrows('nearest', [initpt(1, 1), initpt(1, 2)]);
      end;
    case 89 %Shift - Y
      checkdlg;
      ax = get(gcbf, 'CurrentAxes');
      initpt = get(ax, 'CurrentPoint');
      if ~isempty(initpt)
         arrows('delnearest', [initpt(1, 1), initpt(1, 2)]);
      end;
    case 127 %delete
      cls;
    case 13 %Enter
      if gcbf == i_figno('dialog')
         checkdlg;
         ru;
      end;
    case 6
      i_movie('rewind');
      refresh;
      i_movie('play', 1);
    otherwise
      checkdlg;
      dokeypress(gcbf)
   end;
   
end

function checkdlg
if gcbf == i_figno('dialog')
   close(gcbf);
end;

function resetpoint(H, oldpoint)
if ishandle(H)
   set(H, 'Pointer', oldpoint);
   set(H, 'WindowButtonDown', 'i_callb(''mdown'')');
   set(H, 'WindowButtonMotionFcn', 'i_callb(''mmove'')');
end;
return

function docancel(ix, iy, oldN0, p1, p2, old_Y)
global g_Y;
i_setparvar(ix, p1, []);
i_setparvar(iy, p2, []);
i_keep(transpose(oldN0));
g_Y = old_Y;

function s2 = clearcode(s1)
s2=strrep(s1,'_','');
s2=strrep(s2,'\','');
s2=strrep(s2,'{','(');
s2=strrep(s2,'}',')');
