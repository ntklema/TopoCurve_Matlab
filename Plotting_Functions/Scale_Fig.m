function Scale_Fig(Obj)

p = inputParser;
p.FunctionName = 'Scale_Fig';
addRequired(p,'Obj',@(x) isa(x,'CurveObj'));


end