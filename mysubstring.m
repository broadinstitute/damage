function [s1] = mysubstring(s0,i,j)
s=s0;
if iscell(s0)
    s=char(s0);
end
if (nargin<3)
    j=1;
end
s1=s(:,i:(i+j-1));
if iscell(s0)
    s1=cellstr(s1);
end
