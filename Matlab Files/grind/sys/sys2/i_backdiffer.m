% i_back2 odefile for backward simulation of difference equations
function Result = i_backdiffer(t, Input)
global g_grind TheInput;
TheInput = Input;
Input = 2 * Input - feval(g_grind.odefile, 1, Input);
[Result feval1] = fminsearch(str2func('i_backerror'), Input, optimset('disp','off'));
if feval1 > sum(sum(TheInput)) / 100
   Result = Input * NaN;
end;








