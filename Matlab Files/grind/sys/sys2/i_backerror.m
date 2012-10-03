% i_back2 odefile for backward simulation of difference equations
function Result = i_backerror(Input)
global g_grind TheInput;
Result = sum(abs(TheInput-feval(g_grind.odefile, 1, Input))); % return error

