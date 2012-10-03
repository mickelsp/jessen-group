function [h,isNew] = i_makefig(aname, addnr)
global g_grind;
mnu_null3 = '&3D    (Shift-3)';
mnu_vector = '&Vector field';
mnu_edit = 'Edit &model'; %no convert
mnu_create = 'Create/load &model'; %no convert
switch aname
 case 'changeax'
   domymenu('menubar','axesprop');
   myFig = gcbf;
   oudXlim = g_grind.xaxis.lim;
   oudYlim = g_grind.yaxis.lim;
   oudZlim = g_grind.zaxis.lim;
   pause(0.2);
   waitfor(0, 'CurrentFigure', myFig);
   ax1 = get(myFig, 'CurrentAxes');
   Xlim = get(ax1, 'Xlim');
   Ylim = get(ax1, 'Ylim');
   axchange = 0;
   lab=get(get(gca,'xlabel'),'string');
   nr = i_getno(lab);
   if ~isempty(nr.no) && (lab ~= g_grind.xaxis.var)
      ax('x', lab, Xlim);
      axchange = 1;
   else
      ax('x', [], Xlim);
   end;
   if addnr == 1
      if axchange  ||  max(oudXlim ~= Xlim)  ||  max(oudYlim ~= Ylim)
         if axchange
            close(gcf);
         end;
         null;
      end;
   else
      Ylim = get(ax1, 'Ylim');
      lab=get(get(gca,'ylabel'),'string');
      nr = i_getno(lab);
      if ~isempty(nr.no) && (lab ~= g_grind.yaxis.var)
         ax('y', lab, Ylim);
         axchange = 1;
      else
         ax('y', [], Ylim);
      end;
      if addnr == 2
         if axchange  ||  max(oudXlim ~= Xlim)  ||  max(oudYlim ~= Ylim)
            if axchange
               close(gcf);
            end;
            null;
         end;
      elseif addnr == 3
         Zlim = get(ax1, 'Zlim');
         lab=get(get(gca,'zlabel'),'string');
         nr = i_getno(lab);
         if ~isempty(nr.no) && (lab ~= g_grind.zaxis.var)
            ax('z', lab, Zlim);
            axchange = 1;
         else
            ax('z', [], Zlim);
         end;
         if axchange || max(oudZlim ~= Zlim) || max(oudXlim ~= Xlim) || ...
               max(oudYlim ~= Ylim)
            if axchange
               close(gcf);
            end;
            null3;
         end;
      end;
   end;
   if ~axchange
      set(myFig,'keypressfcn','');
      plotedit('off');
   end;
 case 'daspect'
   das = daspect('mode');
   if strcmp(das, 'manual')
      daspect('auto');
   else
      daspect([1 1 1]);
   end;
 case 'box'
   box;
 case 'grind'
   if ~isempty(g_grind.xaxis)
      iX = i_getno(g_grind.xaxis.var);
      iY = i_getno(g_grind.yaxis.var);
      iZ =  i_getno(g_grind.zaxis.var);
   else
      iX.no = []; iY.no = [];
   end;
   if isempty(g_grind.model)
      set(findobj(gcbf,'tag','mod'),'enable','off');
      set(findobj(gcbf,'tag','mnumodel'),'label',mnu_create);
   else
      set(findobj(gcbf,'tag','mod'),'enable','on');
      set(findobj(gcbf,'tag','mnumodel'),'label',mnu_edit);
   end;
   if (isempty(iX.no) || isempty(iY.no))
      set(findobj(gcbf,'label',mnu_null3),'enable','off')
      set(findobj(gcbf,'label',mnu_vector),'enable','off')
   else
      set(findobj(gcbf,'label',mnu_vector),'enable','on');
      if isempty(iZ.no)
         set(findobj(gcbf,'label',mnu_null3),'enable','off')
      end;
   end;
   if isempty(iY.no) && ~isempty(iX.no)
      if isempty(g_grind) || g_grind.solver.isdiffer
         set(findobj(gcbf, 'tag','mnu2d'),'label','&Iteration map   (Shift-2)','enable','on')
      else
         set(findobj(gcbf, 'tag','mnu2d'),'label','&Differential equation    (Shift-2)','enable','on')
      end;
   elseif isempty(iX.no)
      set(findobj(gcbf, 'tag','mnu2d'),'enable','off')
   else
      set(findobj(gcbf, 'tag','mnu2d'),'label','&2D    (Shift-2)','enable','on')
   end;
 case 'update'
   docheckbox(gcbf);
   docheckhold;
   docheckdasp;
   domymenu('menubar','updatetools',gcbf);
 case 'hold'
   hold;
 otherwise
   if nargin == 1
      addnr = 0;
   end;
   h = figure(i_figno(aname) + addnr);
   isNew = isempty(get(h, 'Children'));
   if isfield(g_grind, 'pen') && ~isempty(g_grind.pen)
      i_plotdefaults(gcf);
   end;
   if nargin == 0
      h = gcf;
   end;
   if isNew  && isfield(g_grind, 'version') 
       %&& (~g_grind.version.isoctave)
      %     hManager = uigetmodemanager(h);
      %     set(hManager.WindowListenerHandles,'Enable','off');
      set(h,'keypressfcn','i_callb(''keypressed'')');
      set(h, 'WindowButtonMotionFcn', 'i_callb(''mmove'')');
      set(h,'ToolBar','figure')
      set(0,'showhiddenhandles','on');
      switch aname
       case 'phase1'
         set(findobj(h,'tag','figMenuToolsAxesProps'),...
            'callback','i_makefig(''changeax'',1)');
       case 'phase2'
         set(findobj(h,'tag','figMenuToolsAxesProps'),...
            'callback','i_makefig(''changeax'',2)');
       case 'phase3'
         set(findobj(h,'tag','figMenuToolsAxesProps'),...
            'callback','i_makefig(''changeax'',3)');
      end;
      lab = '&View mode';
      h2s = findobj(h, 'label', lab);
      if isempty(h2s)
         hgrind=uimenu(h,'label','&GRIND',...
            'callback','i_makefig(''grind'')',...
            'tag','GRIND menu',...
            'pos', 4);
         hru=uimenu(hgrind,'label','&Run',...
            'tag', 'mod','pos', 1);
         uimenu(hru,'label','&Forward (Ctrl-R)',...
            'callback','ru',...
            'pos', 1);
         uimenu(hru,'label','&Backwards  (Ctrl-B)',...
            'callback','backw',...
            'pos', 2);
         uimenu(hgrind,'label','Keep last final state  (Ctrl-K)',...
            'tag', 'mod','callback','ke',...
            'pos', 2);
         uimenu(hgrind,'label','&Time plot  (Ctrl-T)',...
            'tag', 'mod','callback','time',...
            'pos', 3);
         hnull=uimenu(hgrind,'label','&Nullclines',...
            'tag', 'mod','pos', 4);
         uimenu(hnull,'label','&2D    (Shift-2)',...
            'callback','null',...
            'tag','mnu2d',...
            'pos', 1);
         uimenu(hnull, 'label', mnu_null3, ...
            'callback','null3',...
            'pos', 2);
         uimenu(hgrind, 'label', mnu_vector, ...
            'tag', 'mod','callback','vector',...
            'pos', 5);
         uimenu(hgrind,'label','&Erase and make 2D phase plane',...
            'tag', 'mod','callback','e2p',...
            'pos', 6);
         uimenu(hgrind,'label','Continue e&quilibrium (Ctrl-U)',...
            'tag', 'mod','callback','conteq(''-r'');figure(i_figno(''phase2''));',...
            'pos', 7);

         uimenu(hgrind,'label','&Parameter analyser',...
            'tag', 'mod','callback','paranal',...
            'pos', 8);
         uimenu(hgrind,'label','&Find equilibrium',...
            'tag', 'mod','callback','findeq',...
            'separator','on',...
            'pos', 9);
         uimenu(hgrind,'label','&Eigenvalues  (Ctrl-E)',...
            'tag', 'mod','callback','eigen',...
            'pos', 10);
         hsettings=uimenu(hgrind,'label','Settings',...
            'tag','mod','separator','on',...
            'pos', 11);
         uimenu(hgrind, 'label', mnu_edit, ...
            'tag','mnumodel','callback','model',...
            'pos', 12);
         uimenu(hsettings,'label','&Simulation time',...
            'tag', 'mod','callback','simtime',...
            'pos', 1);
         uimenu(hsettings,'label','Set &axes of phase plane',... %No convert
         'tag', 'mod','callback','ax',...
            'pos', 2);
         uimenu(hsettings,'label','&Color of trajectories',...
            'tag', 'mod','callback','setpen',...
            'pos', 3);
         uimenu(hsettings,'label','Set &variables of time plots',... %No convert
         'tag', 'mod','callback','out',...
            'pos', 4);
         uimenu(hsettings,'label','Change parameter',...
            'tag', 'mod','callback','par ?',...
            'pos', 5);
         uimenu(hgrind,'label','&Help GRIND (HTML)',...
            'callback','commands','separator','on',...
            'pos', 12);
         hhelp = findobj(h,'label','&Help');
         if ~isempty(hhelp)
            uimenu(hhelp,'label','&Help GRIND (HTML)',...
               'callback','commands',...
               'pos', 6);
         end;
         htools = findobj(h,'label','&Tools');
         set(htools,'callback','i_makefig(''update'')');
         uimenu(htools,'label','&Box',...
            'callback','i_makefig(''box'')',...
            'pos', 9);
         uimenu(htools,'label','&Hold',...
            'callback','i_makefig(''hold'')',...
            'pos', 10);
         uimenu(htools,'label','&Fixed data aspect ratio',...
            'callback','i_makefig(''daspect'')',...
            'pos', 11);
         uimenu(htools,'label','&Delete last series (Del)',...
            'callback','cls',...
            'pos', 12);
         h2d = uimenu(htools, 'label', lab, ...
            'pos', 2);
         uimenu(h2d,'label','&X vs. Y',...
            'call','set(gca,''View'',[0,90]);',...
            'pos', 1);
         uimenu(h2d,'label','X vs. &Z',... %No convert
         'call','set(gca,''View'',[0,0]);',...
            'pos', 2);
         uimenu(h2d,'label','&Y vs. Z',...
            'call','set(gca,''View'',[90,0]);',...
            'pos', 3);
         uimenu(h2d,'label','&Standard 3D view',...
            'call','set(gca,''View'',[322.5 ,30]);',...
            'pos', 4);
         %     set(hManager.WindowListenerHandles,'Enable','on');
      end;
   end;
   set(0,'showhiddenhandles','off');
end;
function docheckbox(h)
i = get(h, 'Children');
if isempty(i)
   hasbox = 'off';
else
   hasbox = get(gca, 'box');
end;
set(0,'showhiddenhandles','on');
h2s=findobj(gcf,'label','&Box');
set(h2s, 'checked', hasbox);
set(0,'showhiddenhandles','off');
return;
function docheckhold
if ishold
   hashold = 'on';
else
   hashold = 'off';
end;
set(0,'showhiddenhandles','on');
h2s=findobj(gcf,'label','&Hold');
set(h2s, 'checked', hashold);
set(0,'showhiddenhandles','off');
return;
function docheckdasp
das = daspect('mode');
if strcmp(das, 'manual')
   hasdasp = 'on';
else
   hasdasp = 'off';
end;
set(0,'showhiddenhandles','on');
h2s=findobj(gcf,'label','&Fixed data aspect ratio');
set(h2s, 'checked', hasdasp);
set(0,'showhiddenhandles','off');
return;

