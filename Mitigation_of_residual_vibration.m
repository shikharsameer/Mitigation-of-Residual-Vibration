J=5;
% e=(ld^-3*((ld-ld^-5).*(1-(2*ld^2+ld^-4-3)/J)^-1-s1))^0.5;
H=0.05; L=0.05; P=100; G=303000; Br=0.143; u=4*pi*10^(-7); p=2434; %p=rho 
Bapp=7.6646;

c1=2*(H/L)^2;
c=2;
c2=c/(8*H*L*L);
zeta=0.2;
% b=(Br*Bapp/(G*u));
T=(G/(p*L*L));
x1= linspace(0,20,10000);
ld=1.2;
b=((((ld^4)-ld)/(1-((2/ld)+(ld^2)-3)/J)-P*(ld^3)/(4*G*L*L))/ld^3)^0.5;

bl=0;

bu=b;
for i=1:10000
    bs=(bl+bu)/2;
    syms l(t)
    % v1=diff(l,2)*(x*l^(-3)*(2*H*l^(3)+L*L))+y*diff(l,1)-z/(l^4)*(diff(l,2))^2 -m*((G*(l-l^2)/((2*l^(-1)+l^(2)-3)/J)-1)+Bapp*Bo)==0;
    v2=(1+c1*(l^3))*diff(l,2)-((1.5/(l))*(diff(l,1))^2)+3*zeta*c1*(l^3)*diff(l,1)+ 6*(((l^4)-l)/(1-((2/l)+(l^2)-3)/J)-P*(l^3)/(4*G*L*L)-(bs^2)*l^3)==0;
    [V1, S1] = odeToVectorField(v2);
    MMMM=matlabFunction(V1,"vars",["t","Y"]);
    solution=ode45(MMMM,[0 20],[1 0]);
    
    % x1=x*T;
    y1=deval(solution,x1,1);
    y2=deval(solution,x1,2);
   
    if round((bs-bl),4)==0
        for j=1:length(y1)
            if y1(j)==max(y1)
                tau=x1(j);
                bs=bs;
            end
        end
        break
    end
    if max(y1)==ld
        break
    elseif max(y1)>ld
        bu=bs;
    else
        bl=bs;
    end
end


syms l(t)
% v1=diff(l,2)*(x*l^(-3)*(2*H*l^(3)+L*L))+y*diff(l,1)-z/(l^4)*(diff(l,2))^2 -m*((G*(l-l^2)/((2*l^(-1)+l^(2)-3)/J)-1)+Bapp*Bo)==0;
v2=(1+c1*(l^3))*diff(l,2)-((1.5/(l))*(diff(l,1))^2)+3*zeta*c1*(l^3)*diff(l,1)+ 6*(((l^4)-l)/(1-((2/l)+(l^2)-3)/J)-P*(l^3)/(4*G*L*L)-(bs^2)*l^3)==0;
[V1, S1] = odeToVectorField(v2);
MMMM=matlabFunction(V1,"vars",["t","Y"]);
solution=ode45(MMMM,[0 tau],[1 0]);
x1 = linspace(0,tau,10000);
% x1=x*T;
y1=deval(solution,x1,1);
y2=deval(solution,x1,2);
plot(x1,y1,'r','LineWidth',1.5);
hold on

syms l(t)
% v1=diff(l,2)*(x*l^(-3)*(2*H*l^(3)+L*L))+y*diff(l,1)-z/(l^4)*(diff(l,2))^2 -m*((G*(l-l^2)/((2*l^(-1)+l^(2)-3)/J)-1)+Bapp*Bo)==0;
v2=(1+c1*(l^3))*diff(l,2)-((1.5/(l))*(diff(l,1))^2)+3*zeta*c1*(l^3)*diff(l,1)+ 6*(((l^4)-l)/(1-((2/l)+(l^2)-3)/J)-P*(l^3)/(4*G*L*L)-(b^2)*l^3)==0;
[V1, S1] = odeToVectorField(v2);
MMMM=matlabFunction(V1,"vars",["t","Y"]);
solution=ode45(MMMM,[tau 20],[ld 0]);
x1 = linspace(tau,20,10000);
% x1=x*T;
y1=deval(solution,x1,1);
y2=deval(solution,x1,2);
plot(x1,y1,'r','LineWidth',1.5);
grid on
hold on

syms l(t)
% v1=diff(l,2)*(x*l^(-3)*(2*H*l^(3)+L*L))+y*diff(l,1)-z/(l^4)*(diff(l,2))^2 -m*((G*(l-l^2)/((2*l^(-1)+l^(2)-3)/J)-1)+Bapp*Bo)==0;
v2=(1+c1*(l^3))*diff(l,2)-((1.5/(l))*(diff(l,1))^2)+3*zeta*c1*(l^3)*diff(l,1)+ 6*(((l^4)-l)/(1-((2/l)+(l^2)-3)/J)-P*(l^3)/(4*G*L*L)-(b^2)*l^3)==0;
[V1, S1] = odeToVectorField(v2);
MMMM=matlabFunction(V1,"vars",["t","Y"]);
solution=ode45(MMMM,[0 20],[1 0]);
x1 = linspace(0,20,10000);
% x1=x*T;
y1=deval(solution,x1,1);
y2=deval(solution,x1,2);
plot(x1,y1,'g--','LineWidth',1.5);
grid on
hold on

%% Graph from 10-20
disp=y1(length(y1));
vel=y2(length(y2));
J=5;
% e=(ld^-3*((ld-ld^-5).*(1-(2*ld^2+ld^-4-3)/J)^-1-s1))^0.5;
H=0.05; L=0.05; P=100; G=303000; Br=0.143; u=4*pi*10^(-7); p=2434; %p=rho 
Bapp=7.6646;

c1=2*(H/L)^2;
c=2;
c2=c/(8*H*L*L);
zeta=0.2;
% b=(Br*Bapp/(G*u));
T=(G/(p*L*L));
x1= linspace(20,40,10000);
ld=1.5;
b=((((ld^4)-ld)/(1-((2/ld)+(ld^2)-3)/J)-P*(ld^3)/(4*G*L*L))/ld^3)^0.5;

bl=0;

bu=b;
for i=1:10000
    bs=(bl+bu)/2;
    syms l(t)
    % v1=diff(l,2)*(x*l^(-3)*(2*H*l^(3)+L*L))+y*diff(l,1)-z/(l^4)*(diff(l,2))^2 -m*((G*(l-l^2)/((2*l^(-1)+l^(2)-3)/J)-1)+Bapp*Bo)==0;
    v2=(1+c1*(l^3))*diff(l,2)-((1.5/(l))*(diff(l,1))^2)+3*zeta*c1*(l^3)*diff(l,1)+ 6*(((l^4)-l)/(1-((2/l)+(l^2)-3)/J)-P*(l^3)/(4*G*L*L)-(bs^2)*l^3)==0;
    [V1, S1] = odeToVectorField(v2);
    MMMM=matlabFunction(V1,"vars",["t","Y"]);
    solution=ode45(MMMM,[20 40],[disp vel]);
    
    % x1=x*T;
    y1=deval(solution,x1,1);
    y2=deval(solution,x1,2);
   
    if round((bs-bl),4)==0
        for j=1:length(y1)
            if y1(j)==max(y1)
                tau1=x1(j);
                bs=bs;
            end
        end
        break
    end
    if max(y1)==ld
        break
    elseif max(y1)>ld
        bu=bs;
    else
        bl=bs;
    end
end


syms l(t)
% v1=diff(l,2)*(x*l^(-3)*(2*H*l^(3)+L*L))+y*diff(l,1)-z/(l^4)*(diff(l,2))^2 -m*((G*(l-l^2)/((2*l^(-1)+l^(2)-3)/J)-1)+Bapp*Bo)==0;
v2=(1+c1*(l^3))*diff(l,2)-((1.5/(l))*(diff(l,1))^2)+3*zeta*c1*(l^3)*diff(l,1)+ 6*(((l^4)-l)/(1-((2/l)+(l^2)-3)/J)-P*(l^3)/(4*G*L*L)-(bs^2)*l^3)==0;
[V1, S1] = odeToVectorField(v2);
MMMM=matlabFunction(V1,"vars",["t","Y"]);
solution=ode45(MMMM,[20 tau1],[disp vel]);
x1 = linspace(20,tau1,10000);
% x1=x*T;
y1=deval(solution,x1,1);
y2=deval(solution,x1,2);
plot(x1,y1,'r','LineWidth',1.5);
hold on

syms l(t)
% v1=diff(l,2)*(x*l^(-3)*(2*H*l^(3)+L*L))+y*diff(l,1)-z/(l^4)*(diff(l,2))^2 -m*((G*(l-l^2)/((2*l^(-1)+l^(2)-3)/J)-1)+Bapp*Bo)==0;
v2=(1+c1*(l^3))*diff(l,2)-((1.5/(l))*(diff(l,1))^2)+3*zeta*c1*(l^3)*diff(l,1)+ 6*(((l^4)-l)/(1-((2/l)+(l^2)-3)/J)-P*(l^3)/(4*G*L*L)-(b^2)*l^3)==0;
[V1, S1] = odeToVectorField(v2);
MMMM=matlabFunction(V1,"vars",["t","Y"]);
solution=ode45(MMMM,[tau1 40],[ld vel]);
x1 = linspace(tau1,40,10000);
% x1=x*T;
y1=deval(solution,x1,1);
y2=deval(solution,x1,2);
plot(x1,y1,'r','LineWidth',1.5);
grid on
hold on

syms l(t)
% v1=diff(l,2)*(x*l^(-3)*(2*H*l^(3)+L*L))+y*diff(l,1)-z/(l^4)*(diff(l,2))^2 -m*((G*(l-l^2)/((2*l^(-1)+l^(2)-3)/J)-1)+Bapp*Bo)==0;
v2=(1+c1*(l^3))*diff(l,2)-((1.5/(l))*(diff(l,1))^2)+3*zeta*c1*(l^3)*diff(l,1)+ 6*(((l^4)-l)/(1-((2/l)+(l^2)-3)/J)-P*(l^3)/(4*G*L*L)-(b^2)*l^3)==0;
[V1, S1] = odeToVectorField(v2);
MMMM=matlabFunction(V1,"vars",["t","Y"]);
solution=ode45(MMMM,[20 40],[disp vel]);
x1 = linspace(20,40,10000);
% x1=x*T;
y1=deval(solution,x1,1);
y2=deval(solution,x1,2);
plot(x1,y1,'--g','LineWidth',1.5);
ylabel('\lambda',FontSize=12);
xlabel('\tau',FontSize=12);

grid on
hold on


