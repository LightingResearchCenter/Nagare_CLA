 function [tf,inds] = ismonotonic(x)
dx = diff(x);
inds = dx<=0;
tf = ~any(inds);
end