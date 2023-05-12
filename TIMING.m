function [Varb,N]=TIMING(t1,t2,t3,r1,r2,r3)
%see analysis method book,12.1.2, P309
%
M=[r1; r2; r3 ];
T=[t1; t2; t3 ];
mm=M\T;
Varb=1/sqrt(mm(1)^2+mm(2)^2+mm(3)^2);
N=mm*Varb;

