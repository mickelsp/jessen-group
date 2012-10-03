function tilman
global g_grind;
if ~isfield(g_grind, 'tilman');
   error('GRIND:tilman:NoCritVal','Tilman critical values not defined');
end;
[H, new] = i_makefig('phase2');
if new
set(H, 'WindowButtonDown', 'i_callb(''mdown'')');
set(H, 'WindowButtonMotionFcn', 'i_callb(''mmove'')');
end;
set(H, 'Name', 'Phase plane');
oldhold = ishold;
hold('on');
set(gca, 'XLim', g_grind.xaxis.lim);
set(gca, 'YLim', g_grind.yaxis.lim);
xlabel(i_disptext(g_grind.xaxis.var));
ylabel(i_disptext(g_grind.yaxis.var));

critvals = evalin('base', g_grind.tilman.critvals);
rcs = evalin('base', g_grind.tilman.rcs);
maxax = [g_grind.xaxis.lim(2), g_grind.yaxis.lim(2)];
plot([critvals(1, 1), critvals(1, 1), maxax(1)], [maxax(2), critvals(2, 1), critvals(2, 1)], 'b');
plot([critvals(1, 2), critvals(1, 2), maxax(1)], [maxax(2), critvals(2, 2), critvals(2, 2)], 'r');
if all(critvals(:, 1) < critvals(:, 2)) || all(critvals(:, 1) > critvals(:, 2))
   eq1 = critvals(:, 1)';
   eq2 = critvals(:, 2)';
else
   eq1 = max(critvals,[],2);
   eq2 = eq1;
end;
ends = (maxax(1) - eq1(1)) .* rcs(1) + eq1(2);
plot([eq1(1), maxax(1)], [eq1(2), ends(1)], 'b:');
ends = (maxax(1) - eq2(1)) .* rcs(2) + eq2(2);
plot([eq2(1), maxax(1)], [eq2(2), ends(1)], 'r:');
if oldhold
   hold('on');
else
   hold('off');
end;

% function [critval, rc] = getcritval(iW, iR, S, alim, npoints)
% global g_grind;
% oudS = evalin('base', S);
% N0 = i_initvar;
% N0(g_grind.tilman.res) = 9999;
% N0(g_grind.tilman.spec) = 0;
% R = alim(1):(alim(2) - alim(1))  / npoints:alim(2);
% dR = zeros(1, size(R, 2));
% dW = dR;
% N0(iW) = 1;
% for i = 1:size(X, 2)
%    N0(iR) = R(i);
%    assignin('base', S, R(i));
%    yy = feval(g_grind.odefile, 1, N0);
%    dR(i) = yy(iR);
%    dW(i) = yy(iW);
% end;
% assignin('base', S, oudS);
