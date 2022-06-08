%%%%%%%%%%%Figure 2(a)
%%%% condition: w < wt 

clear all
clc

%%%%time initial
t=-0.5e-15:1e-19:5.3333e-15;
M=length(t);

%%%%constant
hb=6.582119514e-16;%unit ev.s

%%%% incident pulse
fc=0.375e15;%for 800 nm
w=2*pi*fc;%carrier wave frequency

%%%%material system
d=0.108;%unit e.nm
a1=0.3; %ratio
a2=0.6;
a3=0.9;
wt1=w/a1;%transition frequency
wt2=w/a2;
wt3=w/a3;

%%%%pulse variable  
E0=0.717997857994181;%unit V/nm
Omegaa=E0*d/hb;
phi=-pi/2;
Cy=32/3*1e-15;%Cycle*2  the pulse length
tao=Cy/4;

%%%%%%computation initilization
h=1e-19;
u=linspace(0,0,M);
v=linspace(0,0,M);
u(1)=1e-12;
v(1)=sqrt(1-u(1)^2);
 
ff1=linspace(0,0,M);
fff1=linspace(0,0,M);
ff2=linspace(0,0,M);
fff2=linspace(0,0,M);
ff3=linspace(0,0,M);
fff3=linspace(0,0,M);
ff1(1)=u(1)/v(1);
fff(1)=abs(ff1(1));
ff2(1)=ff1(1);
fff2(1)=fff1(1);
ff3(1)=ff1(1);
fff3(1)=fff1(1);
 
c1=linspace(0,0,M);
c2=linspace(0,0,M);
cc1=linspace(0,0,M);
cc2=linspace(0,0,M);
ccc1=linspace(0,0,M);
ccc2=linspace(0,0,M); 
c1(1)=u(1);
c2(1)=fff(1)/sqrt(1+fff(1)^2);
 
cc1(1)=c1(1);
cc2(1)=c2(1);
 
ccc1(1)=c1(1);
ccc2(1)=c2(1);
 

 %%%%main loop
    
 %%%% a1=0.3
 %%%% analytical solution
       for k=1:M-1
            theta1=Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*exp(i*wt1*t(k))/(w^2-wt1^2)*(w*sin(w*t(k)+phi)+i*wt1*cos(w*t(k)+phi))-...
            Omegaa*1/2*(square(2*pi/Cy*t(k))+1)/(w^2-wt1^2)*(w*sin(phi)+i*wt1*cos(phi));
            ff1(k+1)=-i/2*(1+exp(-i*(Omegaa*1/2*(square(2*pi/Cy*t(k))+1))^2*wt1/(wt1^2-w^2)*t(k)))*theta1;
            fff1(k+1)=abs(ff1(k+1));
            c2(k+1)=fff1(k+1)/sqrt(1+fff1(k+1)^2);
       end         
 %%%% numerical solution      
        for k=1:M-1
            u1=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(1i*wt1*t(k))*v(k);
            v1=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(-1i*wt1*t(k))*u(k);
            u2=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(1i*wt1*t(k))*(v(k)+0.5*h*v1);
            v2=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(-1i*wt1*t(k))*(u(k)+0.5*h*u1);
            u3=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(1i*wt1*t(k))*(v(k)+0.5*h*v2);
            v3=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(-1i*wt1*t(k))*(u(k)+0.5*h*u2);
            u4=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(1i*wt1*t(k))*(v(k)+h*v3);
            v4=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(-1i*wt1*t(k))*(u(k)+h*u3);
            u(k+1)=u(k)+h*(u1+2*u2+2*u3+u4)/6;
            v(k+1)=v(k)+h*(v1+2*v2+2*v3+v4)/6;
            c1(k+1)=abs(u(k+1));
        end
        
        
          
%%%% a2=0.6
%%%% analytical solution
u=linspace(0,0,M);
v=linspace(0,0,M);
u(1)=1e-12;
v(1)=sqrt(1-u(1)^2);
        for k=1:M-1
            theta2=Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*exp(i*wt2*t(k))/(w^2-wt2^2)*(w*sin(w*t(k)+phi)+i*wt2*cos(w*t(k)+phi))-...
            Omegaa*1/2*(square(2*pi/Cy*t(k))+1)/(w^2-wt2^2)*(w*sin(phi)+i*wt2*cos(phi));
            ff2(k+1)=-i/2*(1+exp(-i*(Omegaa*1/2*(square(2*pi/Cy*t(k))+1))^2*wt2/(wt2^2-w^2)*t(k)))*theta2;
            fff2(k+1)=abs(ff2(k+1));
            cc2(k+1)=fff2(k+1)/sqrt(1+fff2(k+1)^2);
        end
%%%% numerical solution
         for k=1:M-1
            u1=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(1i*wt2*t(k))*v(k);
            v1=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(-1i*wt2*t(k))*u(k);
            u2=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(1i*wt2*t(k))*(v(k)+0.5*h*v1);
            v2=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(-1i*wt2*t(k))*(u(k)+0.5*h*u1);
            u3=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(1i*wt2*t(k))*(v(k)+0.5*h*v2);
            v3=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(-1i*wt2*t(k))*(u(k)+0.5*h*u2);
            u4=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(1i*wt2*t(k))*(v(k)+h*v3);
            v4=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(-1i*wt2*t(k))*(u(k)+h*u3);
            u(k+1)=u(k)+h*(u1+2*u2+2*u3+u4)/6;
            v(k+1)=v(k)+h*(v1+2*v2+2*v3+v4)/6;
            cc1(k+1)=abs(u(k+1));
         end
         
%%%% a3=0.9
%%%% analytical solution
u=linspace(0,0,M);
v=linspace(0,0,M);
u(1)=1e-12;
v(1)=sqrt(1-u(1)^2);
         for k=1:M-1
            theta3=Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*exp(i*wt3*t(k))/(w^2-wt3^2)*(w*sin(w*t(k)+phi)+i*wt3*cos(w*t(k)+phi))-...
            Omegaa*1/2*(square(2*pi/Cy*t(k))+1)/(w^2-wt3^2)*(w*sin(phi)+i*wt3*cos(phi));
            ff3(k+1)=-i/2*(1+exp(-i*(Omegaa*1/2*(square(2*pi/Cy*t(k))+1))^2*wt3/(wt3^2-w^2)*t(k)))*theta3;
            fff3(k+1)=abs(ff3(k+1));
            ccc2(k+1)=fff3(k+1)/sqrt(1+fff3(k+1)^2);
         end       
%%%% numerical solution       
         for k=1:M-1
            u1=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(1i*wt3*t(k))*v(k);
            v1=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(-1i*wt3*t(k))*u(k);
            u2=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(1i*wt3*t(k))*(v(k)+0.5*h*v1);
            v2=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(-1i*wt3*t(k))*(u(k)+0.5*h*u1);
            u3=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(1i*wt3*t(k))*(v(k)+0.5*h*v2);
            v3=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(-1i*wt3*t(k))*(u(k)+0.5*h*u2);
            u4=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(1i*wt3*t(k))*(v(k)+h*v3);
            v4=-1i*Omegaa*1/2*(square(2*pi/Cy*t(k))+1)*cos(w*t(k)+phi)*exp(-1i*wt3*t(k))*(u(k)+h*u3);
            u(k+1)=u(k)+h*(u1+2*u2+2*u3+u4)/6;
            v(k+1)=v(k)+h*(v1+2*v2+2*v3+v4)/6;
            ccc1(k+1)=abs(u(k+1));
         end

%%% transformation
c11=c1./sqrt(1-c1.^2);
c22=cc1./sqrt(1-cc1.^2);
c33=ccc1./sqrt(1-ccc1.^2);
tt=w.*t/(2*pi);

plot(tt,fff1,tt,c11,tt,fff2,tt,c22,tt,fff3,tt,c33)
xlabel('\it{\omegat/2\pi}')
ylabel('| {\itf} |')
