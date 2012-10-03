function i_plotdefaults(H, apen)
global g_grind;
if nargin < 2
   if ~isempty(g_grind) && isfield(g_grind, 'pen')
      apen = g_grind.pen;
   else
      apen = nextpen([]);
   end;
end;
if nargin < 1
   H = gcf;
end;
H=findobj(H,'type','axes');
if isempty(H)
   H = gca;
end;
Cbar=findobj(gcf,'tag','Colorbar');
if ~isempty(Cbar)
   set(Cbar, 'FontSize', apen.fontsize);
end;
set(H, 'FontSize', apen.fontsize);
set(H, 'LineWidth', apen.linewidth);
set(H, 'TickDir', 'out');
%set(H, 'TickDir', apen.tickdir);
set(H, 'Box', apen.box);
t=findobj(H,'type','text');
if ~isempty(t)
   if ~isfield(apen, 'fontsize')
      apen.fontsize = 14;
   end;
   set(t, 'fontsize', apen.fontsize);
end;
if isempty(Cbar) || (gca ~= Cbar)
   View = get(gca, 'view');
   if (apen.fontsize >= 14)  && min(View == [0 90])
      set(gca,'units','normalized');
      if ~isempty(Cbar)
         P = [0.1500    0.1900    0.6455    0.7350];
      else
         P = [0.1300    0.1100    0.7750    0.8150];
      end;
      set(gca, 'position', transform(P,apen.fontsize));
      if ~isempty(Cbar)
         P = [0.831375 0.13 0.058125 0.795];
         set(Cbar, 'position', transform(P,apen.fontsize))
      end;
   end;
end;
function P = transform(P, fontsize)
if fontsize <= 18
   P(2) = P(2) + 0.02;
   P(4) = P(4) - 0.02;
elseif fontsize <= 22
   P(2) = P(2) + 0.04;
   P(4) = P(4) - 0.04;
elseif fontsize <= 28
   P(1) = P(1) + 0.02;
   P(3) = P(3) - 0.02;
   P(2) = P(2) + 0.08;
   P(4) = P(4) - 0.08;
else
   P(1) = P(1) + 0.08;
   P(3) = P(3) - 0.08;
   P(2) = P(2) + 0.16;
   P(4) = P(4) - 0.16;
end;
return

