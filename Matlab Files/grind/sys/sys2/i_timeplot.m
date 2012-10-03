% i_timeplot initial function to create a time plot
function [] = i_timeplot
global g_t g_Y g_grind g_data;
for No = 1:size(g_grind.timevars, 2)
   if ~isempty(g_grind.timevars{No})
      H = i_makefig('time', No - 1);
      plotedit('off');
      oldhold = ishold;
      apen = g_grind.pen;
      apen.i = 1; %  do not cycle colors between runs
      apen.nextpen = 1;
      if ~oldhold
         h=get(get(0,'currentfigure'),'currentaxes');
         if ~isempty(h)
            delete(h);
         end;
      else
         h=get(get(0,'currentfigure'),'currentaxes');
         series = get(h, 'children');
         for i = 1:length(series)
            apen = nextpen(apen);
         end;
      end;
      set(gca,'YLimMode','auto');
      hold('on');
      for i = 1:size(g_grind.timevars{No}, 2)
         ivar = i_varno(g_grind.timevars{No}{i});
         if ~isempty(ivar)
            if ~strcmp(g_grind.outt{No}, 't')
               plot(i_getoutfun(g_grind.outt{No}),g_Y(:, ivar), apen.pen, 'Color', apen.color,'Linewidth',apen.linewidth);
            else
               plot(g_t, g_Y(:, ivar), apen.pen, 'Color', apen.color,'Linewidth',apen.linewidth);
            end;
         else
            if (~isempty(g_data)) && strncmpi(g_grind.timevars{No}{i}, 'Observed ',9)
               ivar2 = i_varno(g_grind.timevars{No}{i}(10:length(g_grind.timevars{No}{i})));
               if ~isempty(g_data) && ~isempty(g_data.obs)
                  plot(g_data.t,g_data.obs(:,ivar2),'o','Color',apen.color,'Linewidth',apen.linewidth);
               else
                  warning('GRIND:time:nodata','No data for "%s"\n', g_grind.timevars{No}{i});
               end;
            else
               res=i_getoutfun(g_grind.timevars{No}{i});
               if ~strcmp(g_grind.outt{No}, 't')
                  tt=i_getoutfun(g_grind.outt{No});
               else
                  tt=g_t;
               end;
               for j=1:size(res,2)   
                  plot(tt,res(:,j), apen.pen,'Color', apen.color,'Linewidth',apen.linewidth);
                  apen=nextpen(apen);
               end;
            end;
         end;
         apen = nextpen(apen);
      end;
      if ~oldhold
         hold('off');
      end;
      tnam = g_grind.diffto;
      f = strfind(tnam, '(');
      tim = tnam(1:f(1) - 1);
      set(H,'Name',[tim 'plot']);
      H1 = get(H, 'CurrentAxes');
      i_plotdefaults(H);
      lim = get(H1, 'Ylim');
      if ~isempty(lim)
         if (lim(2) ~= 0) && ((lim(2) - lim(1)) / lim(2) < 0.001)
            n = 0.01 * lim(2);
            if (lim(1) > n)
               lim(1) = lim(1) - n;
            elseif lim(1) > 0
               lim(1) = 0;
            end;
            lim(2) = lim(2) + n;
            if lim(1)==lim(2)
               lim(2)=lim(1)+eps;
            end;
            set(H1, 'Ylim', lim);
         end
      end;
      if ~strcmp(g_grind.outt{No}, 't')
         %         xlabel([tim '(' g_grind.outt{No} ')']);
         xlabel( g_grind.outt{No});
      else
         xlabel(tnam);
      end;
      nser=length(get(get(0,'currentfigure'), 'children'));
      if nser > 1
         %  if (nser==2)&(
         s = cell(size(g_grind.timevars{No}, 2), 1);
         for i = 1:size(g_grind.timevars{No}, 2)
            s{i} = i_disptext(char(g_grind.timevars{No}{i}));
         end;
         s = i_addlegend(s);
         if length(s) > 1
            legend(s);
            ylabel('');
         else
            if g_grind.version.matlabrelease ~= 11
            legend(char(s{1}));
            legend('hide');
            end;
            ylabel(char(s{1}));
         end;
      elseif size(g_grind.timevars, 2) > 0
         s = i_disptext(g_grind.timevars{No}{1});
         ylabel(s);
         if g_grind.version.matlabrelease ~= 11
            legend(char(s));
            legend('hide');
         end;
      end;
   %   if ~g_grind.version.isoctave
      plotedit('on');
   %   end;
   end;
end;
