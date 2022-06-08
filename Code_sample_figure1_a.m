%%%%%%%%%%%Figure 1(a)

clear all;
clc;

%%%%time initial
t=-20e-15:1e-19:100e-15;
N=length(t);

%%%%constant
hb=6.582119514e-16;% unit eV.s
M=50;

%%%%material system
wba=2.1634e15;%GaAs band gap energy:1.424 eV.
d=0.5; % unit e.nm
gama2=0;

%%%%pulse variable   
E01=0.3454; %unit V/nm 
E02=2.585; 
OmeR1=E01*d/hb;
OmeR2=E02*d/hb;
fai1=1*pi;
%%%%% 5 different phases-second
fai21=-0.5*pi;
fai22=-0.25*pi;
fai23=0*pi;
fai24=0.25*pi;
fai25=0.5*pi;

%%%% incident pulse
tr=50e-15; % time delay
fc=0.375e15;%for 800 nm
wc=2*pi*fc;
env=exp(-(t/(5e-15/1.177)).^2);
env2=exp(-((t-tr)/(5e-15/1.177)).^2);
Ome11=OmeR1.*env.*cos(wc.*t+fai1)+OmeR2.*env2.*cos(wc.*(t-tr)+fai21);
Ome22=OmeR1.*env.*cos(wc.*t+fai1)+OmeR2.*env2.*cos(wc.*(t-tr)+fai22);
Ome33=OmeR1.*env.*cos(wc.*t+fai1)+OmeR2.*env2.*cos(wc.*(t-tr)+fai23);
Ome44=OmeR1.*env.*cos(wc.*t+fai1)+OmeR2.*env2.*cos(wc.*(t-tr)+fai24);
Ome55=OmeR1.*env.*cos(wc.*t+fai1)+OmeR2.*env2.*cos(wc.*(t-tr)+fai25);

% %%%%%%computation initilization
h=1e-19;
u=linspace(0,0,N); %u v w are s1 s2 s3 respectively
v=linspace(0,0,N);
w11=linspace(0,0,N);
w22=linspace(0,0,N);
w33=linspace(0,0,N);
w44=linspace(0,0,N);
w55=linspace(0,0,N);

u(1)=0;
v(1)=0;
w11(1)=-1;

%%%%%%%%main computation_phase 1
for i=1:N-1
    u1=-wba*v(i)-gama2*u(i);
    v1=wba*u(i)+2*Ome11(i)*w11(i)-gama2*v(i);
    w1=-2*Ome11(i)*v(i);
    
    u2=-wba*(v(i)+0.5*h*v1)-gama2*(u(i)+0.5*h*u1);
    v2=wba*(u(i)+0.5*h*u1)+2*Ome11(i)*(w11(i)+0.5*h*w1)-gama2*(v(i)+0.5*h*v1);
    w2=-2*Ome11(i)*(v(i)+0.5*h*v1);
    
    u3=-wba*(v(i)+0.5*h*v2)-gama2*(u(i)+0.5*h*u2);
    v3=wba*(u(i)+0.5*h*u2)+2*Ome11(i)*(w11(i)+0.5*h*w2)-gama2*(v(i)+0.5*h*v2);
    w3=-2*Ome11(i)*(v(i)+0.5*h*v2);
    
    u4=-wba*(v(i)+h*v3)-gama2*(u(i)+h*u3);
    v4=wba*(u(i)+h*u3)+2*Ome11(i)*(w11(i)+h*w3)-gama2*(v(i)+h*v3);
    w4=-2*Ome11(i)*(v(i)+h*v3);
    
    u(i+1)=u(i)+h*(u1+2*u2+2*u3+u4)/6;
    v(i+1)=v(i)+h*(v1+2*v2+2*v3+v4)/6;
    w11(i+1)=w11(i)+h*(w1+2*w2+2*w3+w4)/6;   

end

u(1)=0;
v(1)=0;
w22(1)=-1;

%%%%%%%%main computation_phase 2
for i=1:N-1
    u1=-wba*v(i)-gama2*u(i);
    v1=wba*u(i)+2*Ome22(i)*w22(i)-gama2*v(i);
    w1=-2*Ome22(i)*v(i);
    
    u2=-wba*(v(i)+0.5*h*v1)-gama2*(u(i)+0.5*h*u1);
    v2=wba*(u(i)+0.5*h*u1)+2*Ome22(i)*(w22(i)+0.5*h*w1)-gama2*(v(i)+0.5*h*v1);
    w2=-2*Ome22(i)*(v(i)+0.5*h*v1);
    
    u3=-wba*(v(i)+0.5*h*v2)-gama2*(u(i)+0.5*h*u2);
    v3=wba*(u(i)+0.5*h*u2)+2*Ome22(i)*(w22(i)+0.5*h*w2)-gama2*(v(i)+0.5*h*v2);
    w3=-2*Ome22(i)*(v(i)+0.5*h*v2);
    
    u4=-wba*(v(i)+h*v3)-gama2*(u(i)+h*u3);
    v4=wba*(u(i)+h*u3)+2*Ome22(i)*(w22(i)+h*w3)-gama2*(v(i)+h*v3);
    w4=-2*Ome22(i)*(v(i)+h*v3);
    
    u(i+1)=u(i)+h*(u1+2*u2+2*u3+u4)/6;
    v(i+1)=v(i)+h*(v1+2*v2+2*v3+v4)/6;
    w22(i+1)=w22(i)+h*(w1+2*w2+2*w3+w4)/6;   

end
u(1)=0;
v(1)=0;
w33(1)=-1;

%%%%%%%%main computation_phase 3
for i=1:N-1
    u1=-wba*v(i)-gama2*u(i);
    v1=wba*u(i)+2*Ome33(i)*w33(i)-gama2*v(i);
    w1=-2*Ome33(i)*v(i);
    
    u2=-wba*(v(i)+0.5*h*v1)-gama2*(u(i)+0.5*h*u1);
    v2=wba*(u(i)+0.5*h*u1)+2*Ome33(i)*(w33(i)+0.5*h*w1)-gama2*(v(i)+0.5*h*v1);
    w2=-2*Ome33(i)*(v(i)+0.5*h*v1);
    
    u3=-wba*(v(i)+0.5*h*v2)-gama2*(u(i)+0.5*h*u2);
    v3=wba*(u(i)+0.5*h*u2)+2*Ome33(i)*(w33(i)+0.5*h*w2)-gama2*(v(i)+0.5*h*v2);
    w3=-2*Ome33(i)*(v(i)+0.5*h*v2);
    
    u4=-wba*(v(i)+h*v3)-gama2*(u(i)+h*u3);
    v4=wba*(u(i)+h*u3)+2*Ome33(i)*(w33(i)+h*w3)-gama2*(v(i)+h*v3);
    w4=-2*Ome33(i)*(v(i)+h*v3);
    
    u(i+1)=u(i)+h*(u1+2*u2+2*u3+u4)/6;
    v(i+1)=v(i)+h*(v1+2*v2+2*v3+v4)/6;
    w33(i+1)=w33(i)+h*(w1+2*w2+2*w3+w4)/6;   

end

u(1)=0;
v(1)=0;
w44(1)=-1;

%%%%%%%%main computation_phase 4
for i=1:N-1
    u1=-wba*v(i)-gama2*u(i);
    v1=wba*u(i)+2*Ome44(i)*w44(i)-gama2*v(i);
    w1=-2*Ome44(i)*v(i);
    
    u2=-wba*(v(i)+0.5*h*v1)-gama2*(u(i)+0.5*h*u1);
    v2=wba*(u(i)+0.5*h*u1)+2*Ome44(i)*(w44(i)+0.5*h*w1)-gama2*(v(i)+0.5*h*v1);
    w2=-2*Ome44(i)*(v(i)+0.5*h*v1);
    
    u3=-wba*(v(i)+0.5*h*v2)-gama2*(u(i)+0.5*h*u2);
    v3=wba*(u(i)+0.5*h*u2)+2*Ome44(i)*(w44(i)+0.5*h*w2)-gama2*(v(i)+0.5*h*v2);
    w3=-2*Ome44(i)*(v(i)+0.5*h*v2);
    
    u4=-wba*(v(i)+h*v3)-gama2*(u(i)+h*u3);
    v4=wba*(u(i)+h*u3)+2*Ome44(i)*(w44(i)+h*w3)-gama2*(v(i)+h*v3);
    w4=-2*Ome44(i)*(v(i)+h*v3);
    
    u(i+1)=u(i)+h*(u1+2*u2+2*u3+u4)/6;
    v(i+1)=v(i)+h*(v1+2*v2+2*v3+v4)/6;
    w44(i+1)=w44(i)+h*(w1+2*w2+2*w3+w4)/6;   

end
u(1)=0;
v(1)=0;
w55(1)=-1;

%%%%%%%%main computation_phase 5
for i=1:N-1
    u1=-wba*v(i)-gama2*u(i);
    v1=wba*u(i)+2*Ome55(i)*w55(i)-gama2*v(i);
    w1=-2*Ome55(i)*v(i);
    
    u2=-wba*(v(i)+0.5*h*v1)-gama2*(u(i)+0.5*h*u1);
    v2=wba*(u(i)+0.5*h*u1)+2*Ome55(i)*(w55(i)+0.5*h*w1)-gama2*(v(i)+0.5*h*v1);
    w2=-2*Ome55(i)*(v(i)+0.5*h*v1);
    
    u3=-wba*(v(i)+0.5*h*v2)-gama2*(u(i)+0.5*h*u2);
    v3=wba*(u(i)+0.5*h*u2)+2*Ome55(i)*(w55(i)+0.5*h*w2)-gama2*(v(i)+0.5*h*v2);
    w3=-2*Ome55(i)*(v(i)+0.5*h*v2);
    
    u4=-wba*(v(i)+h*v3)-gama2*(u(i)+h*u3);
    v4=wba*(u(i)+h*u3)+2*Ome55(i)*(w55(i)+h*w3)-gama2*(v(i)+h*v3);
    w4=-2*Ome55(i)*(v(i)+h*v3);
    
    u(i+1)=u(i)+h*(u1+2*u2+2*u3+u4)/6;
    v(i+1)=v(i)+h*(v1+2*v2+2*v3+v4)/6;
    w55(i+1)=w55(i)+h*(w1+2*w2+2*w3+w4)/6;   

end

plot(t*1e15,w11,t*1e15,w22,t*1e15,w33,t*1e15,w44,t*1e15,w55)
xlabel('t/fs')
ylabel('w')