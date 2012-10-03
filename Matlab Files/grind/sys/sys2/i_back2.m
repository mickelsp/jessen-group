% i_back2 odefile for backward simulation of difference equations
function Result = i_back2(t, Input)
global g_grind;
Result = 2 * Input - feval(g_grind.odefile, t, Input);

