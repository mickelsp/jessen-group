function i_makeavi(flag, par1)
global g_frames g_grind t;
switch flag
 case 'addmenu'
    h = findobj(gcf, 'tag', 'GRIND menu');
    if isempty(h)
       h = findobj(gcf, 'tag', 'figMenuTools');
    end;
    if ~isempty(h)
       
   hmenu=uimenu(h,'label','Movie Settings',...
      'tag','viewcells','separator','on',...
      'pos', par1);
   uimenu(hmenu, 'label','Play forward', ...
      'tag','mnumovie','callback','i_movie(''play'',1)',...
      'pos', 1);
   uimenu(hmenu, 'label','Stop', ...
      'tag','mnumovie','callback','i_movie(''deactivate'')',...
      'pos', 2);
   uimenu(hmenu, 'label','Play backward', ...
      'tag','mnumovie','callback','i_movie(''play'',-1)',...
      'pos', 3);
   uimenu(hmenu, 'label','Reset', ...
      'tag','mnumovie','callback','i_movie(''rewind'')',...
      'pos', 4);
   uimenu(hmenu, 'label','Hide buttons', ...
      'tag','hidebuttons','callback','i_movie(''hide'',''off'')',...
      'separator','on','pos', 5);
   if (getrelease~=11) || (exist('CMakeAVI.dll','file') == 3)
      uimenu(hmenu, 'label','Record AVI file', ...
         'tag','avifile','callback','i_movie(''writeavi'')',...
         'pos', 6);
   end;
   end;
 case 'stop/save'
   if isempty(gcbf)
      h = gcf;
   else
      h = gcbf;
   end;
   ud = get(h, 'userdata');
   ud.savingAVI = 0;
   set(h, 'userdata',ud);
   msgbox(['Saved AVI recording to ' ud.filename]);
   if (getrelease == 11)
      aviwrite([ud.pathname ud.filename], g_frames.frames, ud.fps, 'menu');
   else
      g_frames = close(g_frames);
   end;
   mnu=findobj(h,'tag','avifile');
   if ~isempty(mnu)
      set(mnu,'tag','avifile', 'label','Record as AVI file','callback','i_movie(''writeavi'')');
   end;
   dos([ud.pathname ud.filename]);
 case 'init'
   if isempty(gcbf)
      h = gcf;
   else
      h = gcbf;
   end;
   ud = get(h, 'userdata');
   ud.savingAVI = 1;
   if (nargin < 2)||isempty(par1)
      [ud.filename,ud.pathname]=uiputfile('*.avi','Save AVI recording as');
      answer=inputdlg({'Enter the number of frames per second','Start time','End time','The axes only? (Yes/No) (>R11 only)'},'Save to AVI file',...
         1, {'15', num2str(t), num2str(t + g_grind.ndays),'Yes'} );
      ud.fps = str2double(answer{1});
      ud.from = round((str2double(answer{2}) - t) / (g_grind.ndays / g_grind.tstep))+1;
      ud.to = round((str2double(answer{3}) - t) / (g_grind.ndays / g_grind.tstep))+1;
      ud.axisonly = strncmpi(answer{4},'Y',1);
   else
      ud.filename = par1.filename;
      ud.pathname = par1.pathname;
      ud.fps = par1.fps;
      ud.from = par1.from;
      ud.to =par1.to;
      ud.axisonly = par1.axisonly;
   end;
   if isempty(ud.fps)
      ud.fps = 10;
   end;
   if getrelease > 11
      g_frames = avifile([ud.pathname ud.filename]);
      g_frames.Fps = ud.fps;
      g_frames.Quality = 100;
   else
      g_frames.frames = {};
      g_frames.i = 1;
   end;
   %if nargin < 2
  %    msgbox('Recording to AVI file, play (part of) the movie and select "Stop recording" to stop the recording');
  % end;
   set(h, 'userdata',ud);
  % mnu=findobj(h,'tag','avifile');
  % if ~isempty(mnu)
  %    set(mnu,'tag','avifile', 'label','Stop recording','callback','i_movie(''stopavi'')');
  % end;
 case 'getframe'
   if isempty(gcbf)
      h = gcf;
   else
      h = gcbf;
   end;
   if getrelease == 11
      g_frames.frames(g_frames.i) = getindexedframe(gca);
      g_frames.i = g_frames.i + 1;
   else
      ud=get(h,'userdata');
      if ud.axisonly
          F=getframe(gca);
      else
          F = getframe(h);
      end;
      g_frames = addframe(g_frames, F);
   end;
end;


