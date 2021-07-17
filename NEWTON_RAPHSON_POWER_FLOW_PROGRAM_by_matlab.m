% THIS IS THE NEWTON-RAPHSON POWER FLOW PROGRAM
clc;
close all;
clear all;
d2r=pi/180;w=100*pi;
% The Ybus matrix is
[ybus]=ybus;
g=real(ybus);b=imag(ybus);
% The given parameters and initial conditions are
p=[0;-0.96;-0.35;-0.16;0.24];
q=[0;-0.62;-0.14;-0.08;-0.35];
mv=[1.05;1;1;1;1.02];
th=[0;0;0;0;0];
del=1;indx=0;

% The Newton-Raphson iterations starts here
while del>1e-4
end
for i=1:5
temp=0;
for k=1:5
temp=temp+mv(i)*mv(k)*(g(i,k)-1i*b(i,k))*exp(1i*(th(i)-th(k))); 
end
% The mismatches
delp=p-pcal';
delq=q-qcal';
% The Jacobian matrix
for i=1:4
ii=i+1;
for k=1:4
kk=k+1;
j11(i,k)=mv(ii)*mv(kk)*(g(ii,kk)*sin(th(ii)-th(kk))-b(ii,kk)*cos(th(ii)-th(kk)));
end
j11(i,i)=-qcal(ii)-b(ii,ii)*mv(ii)^2;
end
for i=1:4
ii=i+1;
for k=1:4
kk=k+1;
j211(i,k)=-mv(ii)*mv(kk)*(g(ii,kk)*cos(th(ii)-th(kk))-b(ii,kk)*sin(th(ii)-th(kk)));
end
j211(i,i)=pcal(ii)-g(ii,ii)*mv(ii)^2;
end
j21=j211(1:3,1:4);
j12=-j211(1:4,1:3);
for i=1:3
j12(i,i)=pcal(i+1)+g(i+1,i+1)*mv(i+1)^2;
end
j22=j11(1:3,1:3);
for i=1:3
j22(i,i)=qcal(i+1)-b(i+1,i+1)*mv(i+1)^2;
end
jacob=[j11 j12;j21 j22]; 
delpq=[delp(2:5);delq(2:4)]; 
corr=inv(jacob)*delpq; 
th=th+[0;corr(1:4)];
mv=mv+[0;mv(2:4).*corr(5:7);0]; 
del=max(abs(delpq));
indx=indx+1;
end
preal=(pcal+[0 0 0 0 0.24])*100;
preac=(qcal+[0 0 0 0 0.11])*100;
% Power flow calculations
for i=1:5
v(i)=mv(i)*exp(1i*th(i));
end
for i=1:4
for k=i+1:5
if (ybus(i,k)==0)
s(i,k)=0;s(k,i)=0;
c(i,k)=0;c(k,i)=0;
q(i,k)=0;q(k,i)=0;
cur(i,k)=0;cur(k,i)=0;
else
cu=-(v(i)-v(k))*ybus(i,k);
s(i,k)=-v(i)*cu'*100;
s(k,i)=v(k)*cu'*100;
c(i,k)=100*abs(ybus(i,k))*abs(v(i))^2;
c(k,i)=100*abs(ybus(k,i))*abs(v(k))^2;
cur(i,k)=cu;cur(k,i)=-cur(i,k);
end
end
end
pwr=real(s);
qwr=imag(s);
q=qwr-c;
% Power loss
ilin=abs(cur);
for i=1:4
for k=i+1:5
if (ybus(i,k)==0)
pl(i,k)=0;pl(k,i)=0;
ql(i,k)=0;ql(k,i)=0;
else
z=-1/ybus(i,k);
r=real(z);
x=imag(z);
pl(i,k)=100*r*ilin(i,k)^2;pl(k,i)=pl(i,k); 
ql(i,k)=100*x*ilin(i,k)^2;ql(k,i)=ql(i,k);
end
end
end