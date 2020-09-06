function ek = nhieu(t)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
ek=0;
for w=1:10
    ek=ek+0.06/w*sin(w*t);
end
end