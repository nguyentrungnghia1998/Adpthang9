Q1d=0.001;
Rd=0.001;
gamma=0.9;
Q2d=100;
siold=[0 0.2 0.6 0.6]';
yold=[0 0.2 0.6]';
rold=[0.4 0.3]';
y_=[0 0.2 0.6 0 0.2 0.6]';
r_=[0.4 0.3 0.4 0.3]';
Pd=ones(12)*0.0001;
ek=nhieu(0);
unew=-1/(Rd/gamma+Pd(1,1))*(Pd(2)*0+Pd(1,3:8)*y_+Pd(1,9:12)*r_)+ek;
y2=mohinh(siold,unew,rold);