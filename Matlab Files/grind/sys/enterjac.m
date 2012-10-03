%ENTERJAC   Enter equations of the Jacobian matrix
%   Enter the analytical Jacobian matrix, consisting of the derivatives
%   of the right-hand sides of all differential/difference equations with respect to each state
%   variable. 
%   This is used to calculate the linearized system in equilibria. This way
%   you can study the stability of equilibria. 
%
%   Usage:
%   ENTERJAC - the user is prompted to give the elements of the Jacobian 
%   ENTERJAC({ J11,J12;J21,J22}) - the whole matrix is entered in one message 
%   (J11 is the first string with an equation. 
%   ENTERJAC 1 1 J11 - each element is entered separarely
%
%   See also eigen, findeq, lyapspect

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function enterjac(varargin)
global g_grind;
i_parcheck;
n = g_grind.statevars.dim;
if isempty(g_grind.Jacobian)
    g_grind.Jacobian = cell(n, n);
end
if nargin == 1
   g_grind.Jacobian=varargin{1};
elseif nargin==3
  i=i_checkstr(varargin{1});
  j=i_checkstr(varargin{2});
  g_grind.Jacobian{i,j}=varargin{3};
elseif nargin == 0
   if n == 0
      error('GRIND:enterjac:NoStatevars','No state variables');
   elseif n>20
      error('GRIND:enterjac:TooManyStatevars','Too many state variables')
   end;
   prompt = cell(1, n);
   def=prompt;
   hasempty = 0;
   for i = 1:n
      g_rhs=i_statevars_names(i);
      for j = 1:n
         prompt{j}=sprintf('Enter partial derivative of [ %s'' ] with respect to %s',g_rhs ,i_statevars_names(j));
         def{j} = char(g_grind.Jacobian{i, j});
      end
      answer = inputdlg(prompt, 'Enter Jacobian', 1, def);
      if ~isempty(answer)
         for j = 1:n
            g_grind.Jacobian{i, j} = answer{j};
            if isempty(answer{j})
               hasempty = hasempty + 1;
            end;
         end
      else
         if isempty(g_grind.Jacobian{1})
            g_grind.Jacobian = [];
         end;
      end
   end
   if hasempty == n * n
      g_grind.Jacobian = [];
      error('GRIND:enterjac:EmptyJac','The entered Jacobian has only empty cells');
   end;
end;
if n == 1
   fprintf('J=[''%s'']\n',g_grind.Jacobian{1,1});
elseif n == 2
   fprintf('J=[%s, %s;\n %s, %s]\n',g_grind.Jacobian{1,1},g_grind.Jacobian{1,2}, ...
      g_grind.Jacobian{2, 1}, g_grind.Jacobian{2, 2});
end;
