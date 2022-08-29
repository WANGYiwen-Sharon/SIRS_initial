%% odef_simple_sirs_protection_nis_par
%% add natural infection (ni)
%% add full susceptible people (s)
%% add secondary/additional boosting (b)
%% Simple SIRS with: 1. evolving virus
%%                   2. antibody boosting with natural infection
%% Can pass parameters
function xdot = odef_simple_sirs_protection_ni_par(t,x,par)
global mu beta gamma ci wan;
%mu: birth and death rate
%beta: infection force
%gamma: recovery rate
%mu = 0;
%beta = 0.497*0.7;
%beta = 0.3636;
%gamma = 1/5;
%N = 1000;

s = x(1); % Suscpetible after waning immunity
i = x(2); % Infected 
v1 = x(3); % Transient immunity
sr = x(4);% Partial protection before waning immunity
a = x(5); % Population antigenicity
v = x(6); % Virus antigenicity
l = x(7); % Deleterious effect

% 1: Naive group
% 2: Vaccinated group
% 3: Non-vaccinated group withous previous infection
% 4: Vaccinated group with previous infection

% s0: Naive group
s0_n = x(8);
i_n = x(9);
r_n = x(10);
sr_n = x(11);
sl_n = x(12);
ci_x = x(13);

% r_n: A longer transient immunity for natural infection
% sr_n: Partial protection before waning immunity for natural infection
% s_n: Susceptible after waning immunity for natural infection


% https://www.hopkinsmedicine.org/health/conditions-and-diseases/coronavirus/covid-natural-immunity-what-you-need-to-know
% https://www.ettoday.net/news/20220120/2173425.htm
% https://www.health.ny.gov/press/releases/2022/2022-01-19_covid_reinfection_study.htm
% For example, in New York State in the week beginning October 3, unvaccinated people without a previous COVID-19 diagnosis, were 4.5 fold more likely to have a positive COVID-19 test than vaccinated people without a previous COVID-19 diagnosis (Group 1 vs. 2), 14.7-fold more likely than unvaccinated people with a previous diagnosis (Group 1 vs. 3), and 19.8-fold more likely than vaccinated people with a previous COVID-19 diagnosis (Group 1 vs. 4).
% Vaccine induced immunity protect 77.8% 
% Infection induced immunity cover 93.2%; 1.2 fold than vaccine-induced
% immunity

d = v-a;
d_vac = 1;
if d_vac > d
  d_vac = d;
end
%l = 0.2;
pr = 1 - d - l; %protection is 1 minus antigenic distance and deleterious effect
if pr >= 1
    pr = 1;
end
%pr = 1 - v ;
alpha = 0;
% vaccination
%Fixed vaccine doses
%50%
v_cov = 0.6; % Default 0.6; > 0.08 or < 0.02 biannual; 0.05 becomes once a year
delay = 0;


%if t > 365+delay && t < 365+30*13+delay
%    alpha = v_cov/(30*13);
%end

% 未來持續注射 Vaccination not stop
if t > 365+delay
    alpha = v_cov/(30*13);
end
% 未來注射減少 Vaccination reduces
if t > 365*2
  alpha = 1*v_cov/(30*13);
end
%delay = delay + 180;
%if t>400+delay && t <=460+delay;
    %alpha = 0.001;
%    alpha = v_cov/60; % vaccinatation daily rate
    %if r > 0.6
    % alpha = 0;
    %end
%end
% S' = mu - beta*S*I - mu*S
% I' = beta*S*I - gamma*I - mu*I;
% R' = gamma*I - mu*R

%s_dot = mu - beta*s.*i - mu*s;
%i_dot = beta*s.*i - gamma*i - mu*i;
%r_dot = gamma*i - mu*r; 
%xdot = [s_dot; i_dot; r_dot];

% Four important factors determine the future transition from pandemic to
% epidemic
% 1. Protectiveness after loss of antibody (Loss)
% 2. Virus mutation (Mut)
% 3. Vaccination coverage (V_cov)
% 4. Short term immunity period (Imm_dur)

Loss = 0.8;     %Low 0.2; Medium 0.5; High 0.8
Loss_ni = Loss*0.95; 
Mut = 1;        %Mutation speed: 1
if t > 365*1.9 && t < 365*2.1
%    Mut = 1; % 1 -> 2
end
Imm_dur = 42;   %Transient immunity: 42 days default
% 自然免疫的Immune duration should be longer
Mob = -0.2;        %Change in mobility
Amp = 0.35;        % 0.35 Seasonal, 0.1 Endemic
if t < 365*1.5
 Mob = -0.4;
end
if t > 365 && 5 < 365*1.8
% Mob = -0.3;
end
if t > 365*1.8
 Mob = 0;
end

if t > 365*10
   pr; 
end
pr_ni = pr*1.2;
if pr_ni > 1
  pr_ni = 1;
end
cont = Amp*sin(1/180*(pi*2*(t+45)))+1+Mob;
%cont should be google mobility x (1 - )

%Transmission rate of naturally infected individuals
beta_naive = beta*cont;
beta_eff_ni = beta*cont*(1-pr_ni); % NI protects 1.2 fold than vaccination
beta_eff_loss_ni  = beta*cont*(1-pr_ni*(1-Loss));

%Transmission rate of vaccinated individuals
beta_eff = beta*cont*(1-pr);
beta_eff_loss = beta*cont*(1-pr*(1-Loss));

%%
i_total = i+i_n;
%s(waning immunity) r(fully recovere) sr(vaccine protected)
%s_dot = -beta_eff_loss*(1-i-r)*i - alpha*(1-i-r) + wan*r;
s_dot = -beta_eff_loss*s*i_total - alpha*s + wan*sr; % s = 1-i-r-sr
i_dot = beta_eff_loss*s*i_total + beta_eff_loss*v1*i_total + beta_eff*sr*i_total - gamma*i; % i -> natural infection  group
%change r_dot to double boosting
sr_dot = alpha*(s + sr + sl_n + sr_n) - beta_eff*sr*i_total - alpha*sr - wan*sr;
v1_dot = alpha*s0_n - wan*v1 - beta_eff_loss*v1*i_total;

% natural infection
s0_n_dot = -beta_naive*s0_n*i_total - alpha*s0_n + wan*v1;
i_n_dot = beta_naive*s0_n*i_total + beta_eff_loss_ni*sl_n*i_total + beta_eff_ni*sr_n*i_total - gamma*i_n;
r_n_dot = gamma*i + gamma*i_n - 1/Imm_dur*r_n;
sr_n_dot = 1/Imm_dur*r_n - beta_eff_ni*sr_n*i_total - alpha*sr_n - wan*sr_n;
sl_n_dot = -beta_eff_loss_ni*sl_n*i_total - alpha*sl_n + wan*sr_n; % s = 1-i-r-sr

% Every new infections increase antibody protection
% The rate of antibody protection is proportional to distance
a_dot = d*beta_eff_loss*s*i_total + d*beta_eff_loss*v1*i_total + d*beta_eff*sr*i_total + d*beta_naive*s0_n*i_total + d*beta_eff_loss_ni*sl_n*i_total + d*beta_eff_ni*sr_n*i_total + d*alpha*(s+sr+sl_n+sr_n); % Only naive population gain immunity
% Virus evolution
v_dot = Mut*((pr*(1-Loss)*beta_eff_loss*s*i_total + pr*(1-Loss)*beta_eff_loss*v1*i_total + pr*beta_eff*sr*i_total) + 0*beta_naive*s0_n*i_total + (pr_ni*(1-Loss)*beta_eff_loss_ni*sl_n*i_total + pr_ni*beta_eff_ni*sr_n*i_total));
l_dot = 0;
%l_dot = 0.00001*(v_dot - v*(beta_eff_loss*(1-i-r)*i + beta_eff*r*i));
ci_dot = beta_eff_loss*v1*i_total + beta_eff_loss*s*i_total + beta_naive*s0_n*i_total + beta_eff_loss_ni*sl_n*i_total + beta_eff*sr*i_total + beta_eff_ni*sr_n*i_total;
xdot = [s_dot; i_dot; v1_dot; sr_dot; a_dot; v_dot; l_dot; s0_n_dot; i_n_dot; r_n_dot; sr_n_dot; sl_n_dot; ci_dot];

end




