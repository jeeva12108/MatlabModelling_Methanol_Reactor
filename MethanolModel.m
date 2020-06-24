%{
Methanol Reactor Modelling
Authors-Jeeva, Gowtham, Balasubramanian
Project at IIT Madras
Adiabatic Plot of Methanol synthesis from Sync Gas
Kinetic Model: Arhenius Equation
}%

function Fun=conversion_all_froment_with_delh(~,Fin)
global P R rho A eb
Ftotmol=sum(Fin(1:6));
T=Fin(7);
xco2=(Fin(1)/Ftotmol);
xh2=(Fin(2)/Ftotmol);
xmeoh=(Fin(3)/Ftotmol);
xh2o=(Fin(4)/Ftotmol);
xco=(Fin(5)/Ftotmol);
xinert=(Fin(6)/Ftotmol);
pco2=xco2*P;
ph2=xh2*P;
pmeoh=xmeoh*P;
ph2o=xh2o*P;
pco=xco*P;
K2=34553.38;
K3=0.499*exp(17197/(R*T));
K4=6.62e-11*exp(124119/(R*T));
E=1+(K2*(ph2o/ph2))+((K3)*(sqrt(ph2)))+(K4*(ph2o));%Arhenius Model 
reactcoeff1=1.07*exp(36696/(R*T));
k1eq=10^((3066/T)-10.592);
r1=(reactcoeff1*pco2*ph2*(1-
((ph2o*pmeoh)/(k1eq*((ph2^3))*pco2))))/((E^3));
delh10=-49510;
k1=1.22e10*exp((-94765)/(R*T));
k3eq=10^((2073/T)-2.029);
r2=((k1*pco2)*(1-(k3eq*ph2o*pco/(pco2*ph2))))/(E);
delh20=41190;
cpco2=R*(5.457+(1.045*(10^(-3)))*T-(1.157*(10^5))*(T^(-
2)));

cph2=R*(3.249+(0.422*(10^(-3)))*T+(0.083*(10^5))*(T^(-
2)));
cpmeoh=R*(2.211+(12.216*(10^(-3)))*T-(3.45*(10^(-
6)))*(T^2));
cph2o=R*(3.470+(1.45*(10^(-3)))*T+(0.121*(10^5))*(T^(-
2)));
cpco=R*(3.376+(0.557*(10^(-3)))*T-(0.031*(10^5))*(T^(-
2)));
cpinert=20.786;
cpt=cpco2*xco2+cph2*xh2+cpmeoh*xmeoh+cph2o*xh2o+cpco*xco+c
pinert*xinert;
dela1=R*(2.211+3.470-(3*3.249)-5.457);
delb1=R*(12.216*(10^(-3))+1.45*(10^(-3))-(3*0.422*(10^(-
3)))-1.045*(10^(-3)));
delc1=-R*(3.45*(10^(-6)));
deld1=R*((0.121*(10^5))+(1.157*(10^5))-3*(0.083*(10^5)));
To=298.15;
delh1=delh10+dela1*(T-To)+(delb1/2)*((T^2)-
(To^2))+(delc1/3)*((T^3)-(To^3))+deld1*((1/To)-(1/T));
dela2=R*(3.376+3.470-3.249-5.457);
delb2=-R*(1.45*(10^(-3))+0.577*(10^(-3))-1.045*(10^(-3))-
0.422*(10^-(3)));
delc2=0;
deld2=-R*((0.031*(10^5))+(0.121*(10^5))+(1.157*(10^5))-
(0.083*(10^5)));
delh2=delh20+dela2*(T-To)+(delb2/2)*((T^2)-
(To^2))+delc2/3*((T^3)-(To^3))+deld2*((1/To)-(1/T));
Fun(1)=-(r1+r2)*rho*A*(1-eb);
Fun(2)=-(r1*3+r2)*rho*A*(1-eb);
Fun(3)=r1*rho*A*(1-eb);
Fun(4)=(r1+r2)*rho*A*(1-eb);
Fun(5)=r2*rho*A*(1-eb);
Fun(6)=0*rho*A*(1-eb);
Fun(7)=rho*(1-eb)*A*(r1*(-1*delh1)+r2*(-
1*delh2))/(Ftotmol*cpt);
Fun=Fun&#39;;

end

clc
clearall
closeall
T=493;
global P R rho A D eb
P=50;
rho=1775;
m_cat=0.0348;
v_cat=m_cat/rho;
D=0.016;
L=1.5;
A=(pi*(D^2))/4;
eb=1-(v_cat/(A*0.15));
R=8.314;
m_flo=2.8*(10^-5);
xco2=0.03;
xh2=0.82;
xmeoh=0;
xh2o=0;
xco=0.04;
xinert=0.11;
mco2=0.04401;
mh2=0.0020158;
mmeoh=0.0324;
mh2o=0.018;
mco=0.02801;
minert=0.03995;
moverall=(mco2*xco2)+(mh2*xh2)+(mmeoh*xmeoh)+(mh2o*xh2)
+(mco*xco)+(minert*xinert);
Ftotmol=m_flo/moverall;
Fco2=xco2*Ftotmol;
Fh2=xh2*Ftotmol;
Fmeoh=xmeoh*Ftotmol;
Fh2o=xh2o*Ftotmol;
Fco=xco*Ftotmol;
Finert=xinert*Ftotmol;
z= 0 : 0.02 : L;
Fin=[Fco2 Fh2 Fmeoh Fh2o FcoFinert T];
[z,Fout]=ode15s(@conversion1,z,Fin);
figure(1); %Ploting
plot(z/L,(Fout(:,1)/Ftotmol)*100)

xlabel(&#39;z/L(m)&#39;)
ylabel(&#39;mole% of co2&#39;)
figure(2);
plot(z/L,(Fout(:,2)/Ftotmol)*100)
xlabel(&#39;z/L(m)&#39;)
ylabel(&#39;mole% of h2&#39;)
figure(3);
plot(z/L,(Fout(:,3)/Ftotmol)*100)
xlabel(&#39;z/L(m)&#39;)
ylabel(&#39;mole% of meoh&#39;)
figure(4);
plot(z/L,(Fout(:,4)/Ftotmol)*100)
xlabel(&#39;z/L(m)&#39;)
ylabel(&#39;mole% of h2o&#39;)
figure(5);
plot(z/L,(Fout(:,5)/Ftotmol)*100)
xlabel(&#39;z/L(m)&#39;)
ylabel(&#39;mole% of co&#39;)
figure(6);
plot(z/L,Fout(:,7))
xlabel(&#39;z/L(m)&#39;)
ylabel(&#39;temp (K)&#39;)
