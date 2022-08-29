function main()

disp('current theta = ');
%disp(theta)
% observed data at time point 1, 2, 3, ..., 6
%data = [10.^2.6 10.^5.0 10.^5.1 10.^4.9 10.^3.8 10.^1.9];

%global k p c d;

global mu beta gamma ci wan A0 V0;
%beta = 0.3636;
%beta = 0.45*2;
gamma = 1/6;
R_0 = 2;
beta = R_0 * gamma;
wan = 1/180; %1/180
%VacCov = 0.67; % vaccine coverage percentage
VacCov = 0;
if ~exist('I0') % set initial infecteds I0
  I0 = 0.00001 
end
if ~exist('S0')
  S0 = 1-I0-VacCov;
end

R0 = 0;
%SR0 = 1-S0-I0-R0;
SR0 = 0;
SL0 = 0;
A0 = 0; % Host Antigenicity
V0 = 1; % Virus Antigenicity
L0 = 0; % Deleterious Effect
CI = 0;
% single strain
simple_sirs();
%pause


function simple_sirs()
starttime = 0;
totaltime = 365*6.5*2;
t = [starttime:1:totaltime];
ci = 0;
par.beta = R_0 * gamma;
par.Mut = 1;
par.Imm_dur = 42; %(days)
par.RepRate = 0.8; %(reporting rate) 
par.I50 = 1000000; 

%[t1, x] = ode45('odef_simple_sirs_protection_c',t,[0; 0; 0; 0; A0; V0; L0; S0; I0; R0; SR0]);
[t1, x] = ode45(@odef_simple_sirs_protection_nisb_par,t,[0; 0; 0; 0; A0; V0; L0; S0; I0; R0; SR0; SL0; CI],[],par);

    
figure;
hold;
plot(t1,x(:,9),'r','LineWidth',1.2);
plot(t1,x(:,2),'b','LineWidth',1.2);
legend({'Infected in previously native individuals','Infected in previously vaccinated individuals'},'FontSize',16);
xticks([0 182.5 182.5*2 182.5*3 182.5*4 182.5*5 182.5*6 182.5*7 182.5*8 182.5*9 182.5*10 182.5*11 182.5*12 182.5*13 182.5*14])
xticklabels({'2020W','2020S','2021W','2021S','2022W','2022S','2023W','2023S','2024W','2024S','6W','6S','7W','7S'})
axisHandle = gca;
axisHandle.XAxis.TickLabelRotation = 90;


figure;
hold;
plot(t1,x(:,4),'r','LineWidth',1.2);
plot(t1,x(:,11),'b','LineWidth',1.2);
plot(t1,x(:,10),'k','LineWidth',1.2);
%Y = [x(:,8) x(:,12) x(:,1)]; %Naive, 
%area(Y,'LineStyle','none');
%newcolors = [0.1 0.6 1; 0.1 0.7 1; 0.1 0.8 1];
%colororder(newcolors)
%%b8 = plot(tt(t0:end),local_total_mean7(t0:end),'color',[101/255 171/255 9/255],'linestyle',':','linewidth',2.5)
%h8=fill(XX,YY,[237/255 218/255 218/255],'Linestyle','none');
legend({'Vaccine protected','Natural infection protected','Full Protection'},'FontSize',16)
xticks([0 182.5 182.5*2 182.5*3 182.5*4 182.5*5 182.5*6 182.5*7 182.5*8 182.5*9 182.5*10 182.5*11 182.5*12 182.5*13 182.5*14])
xticklabels({'2020W','2020S','2021W','2021S','2022W','2022S','2023W','2023S','2024W','2024S','6W','6S','7W','7S'})
axisHandle = gca;
axisHandle.XAxis.TickLabelRotation = 90;


figure;
hold;
plot(t1,x(:,8),'r');
plot(t1,x(:,12),'b');
plot(t1,x(:,1),'k');
legend({'Susceptible','Naturally infected after waning','Vaccinated after waning'},'FontSize',16)
xticks([0 182.5 182.5*2 182.5*3 182.5*4 182.5*5 182.5*6 182.5*7 182.5*8 182.5*9 182.5*10 182.5*11 182.5*12 182.5*13 182.5*14])
xticklabels({'2020W','2020S','2021W','2021S','2022W','2022S','2023W','2023S','2024W','2024S','6W','6S','7W','7S'})
axisHandle = gca;
axisHandle.XAxis.TickLabelRotation = 90;


%plot(t1,1-x(:,1)-x(:,2),'k');
mortality_ratio = x(:,3)./(x(:,2)+x(:,3));
diff_mr50 = abs(mortality_ratio - 0.7);
date_mr50 = find(diff_mr50 == min(diff_mr50));
date_mr50
infected_ratio = x(:,2)+x(:,3);
infected_ratio(end,1)

A = x(:,5);
V = x(:,6);
L = x(:,7);
D = V - A;
Prot = 1 - D;
%plot(t1,V-A, 'b');
xticks([0 182.5 182.5*2 182.5*3 182.5*4 182.5*5 182.5*6 182.5*7 182.5*8 182.5*9 182.5*10 182.5*11 182.5*12 182.5*13 182.5*14])
xticklabels({'2020W','2020S','2021W','2021S','2022W','2022S','2023W','2023S','2024W','2024S','6W','6S','7W','7S'})
axisHandle = gca;
axisHandle.XAxis.TickLabelRotation = 90;


Loss = 0.8; % check ode function
Avg_prot = x(:,4).*Prot + x(:,1).*Prot.*(1-Loss);
%plot(t1, Avg_prot);
%plot(t1,D);
figure;
hold;
plot(t1,A, 'r');
plot(t1,V, 'b');

legend({'people antigenicity','virus antigenicity'},'FontSize',16)

%plot(t1,V);
%hold;
%plot(t1,A);
leng = length(x(:,13))
ci_last = x(leng,13);
ci_2yr = x(floor(365*2.2),13);
ci_annual = (ci_last-ci_2yr)/((leng-(365*2.2))/365)
ci_2022w = (x(floor(365*2.2),13)-x(floor(365*1.8),13))

%每月平均
clear all;

end


function simple_simi()
t = [0:1:100];
[t1, x] = ode45('odef_simple_simi',t,[1; 0.01; 0; 0]);
figure;
plot(t1,x(:,1));
hold;
plot(t1,x(:,2),'r');
plot(t1,x(:,3),'g');
plot(t1,x(:,4),'c');
hleg1 = legend('S','I','M','I2');
end


function two_strain_simi()
t = [0:1:200];
[t1, x] = ode45('odef_twostrain_simi',t,[1; 0.01; 0; 0; 0; 0.005; 0; 0; 0; 0; 0; 0]);

Nss = x(:,1);
Nis = x(:,2);
Nms = x(:,3);
Ni2s = x(:,4);
Nmi = x(:,5);
Nsi = x(:,6);
Nsm = x(:,7);
Nsi2 = x(:,8);
Nim = x(:,9);
Nmm = x(:,10);
Nmi2 = x(:,11);
Ni2m = x(:,12);


Nig = Nis+Nim;
Nmg = Nms+Nmi+Nmi2+Nmm;
Ni2g = Ni2s+Ni2m;

figure;
plot(t1,x(:,1));
hold;
plot(t1,Nig,'r');
plot(t1,Nmg,'g');
plot(t1,Ni2g,'k');
%plot(t1,Nmm,'y');
hleg1 = legend('Nss','Nig','Nmg','Ni2g');
end


end
