# SIR Protection nisb
## 1. Model Structure
![alt text](./model%20structure.png)
## 2. main_sir_protection_nisb Global Parameters
| Parameters | Value | Comments |
| :--------: | :---: | :-----: |
| gamma  | 1/6 | removal rate |
| R_0  | 2 | |
| beta | R_0 * gamma | transmission rate |
| wan | 180 | antibody vanishing rate |
| VacVov | 0 | initial vaccine coverage |
| I0 | 0.00001 | initial infected |
| S0 | 1 - I0 - VacCov | initial susceptible |
| R0 | 0 | recovered | 
| SR0 | 0 | 1 - S0 - I0 - R0 initial antibody protected population |
| SL0 | 0 | initial antibody vanishing population |
| A0 | 0 | host antigencity |
| V0 | 1 | various antigencity |
| L0 | 0 | deleterious effect |
| CI | 0 | get infected population |
| simple sirs() |
| starttime | 0 | start time of the model |
| totaltime | 365 * 6.5 * 2 | total time period of the model |
| t | [starttime:1:totaltime] | time period |
| ci | 0 | |
| par.beta | R0 * gamma | transmission rate |
| par.Mut | 1 | mutation speed |
| par.Imm_dur | 42 | short term immunity period ($R \Rightarrow S_{R}$ )| 
| par.RepDate | 0.8 | reporting rate |
| par.I50 | 1000000 | |

## 3. simple_sirs_protection_nisb_par Global Parameters
| Parameters | Value | Comments |
| :--------: | :---: | :-----: |
| s | x(1) | susceptible after waning immunity ($V_{L}$)|
| i | x(2) | infected from &V_1, V_L, V_F&($I$)|
| v1 | x(3) | first dose of vaccine ($V_1$) |
| sr | x(4) | partial protection before waning immunity ($V_F$)
| a | x(5) | population antigencity |
| v | x(6) | virus antigencity |
| l | x(7) | deleterious effect | 
| S0_n | x(8) | ordinary suspectible without any antibody ($S_0$)
| i_n | x(9) | infected from $S_0, S_L, S_N$ ($I_n$) |
| r_n | x(10) | recovered population with full antibody ($R$) |
| sr_n | x(11) | recovered population with first weakening of antibody($S_R$) |
| sl_n | x(12) | recovered population with second weakening of antibody ($S_L$)|
| ci_c | x(13) | get infected population |
| |
| d | v - a | antigenitic distance | 
| d_vac | min(d, d_vac) | |
| alpha | v_cov/(30*13) when t > 365 - delay <br> 1 * v_cov/(30*13) when t > 365 * 2 | vaccine doses rate |
| v_cov | 0.6 | vaccine coverage |
| delay | 0 | vaccination delay |
| Loss | 0.8 | Protectiveness after loss of antibody <br> (low 0.2 medium 0.5 high 0.8) |
| Loss_ni | Loss * 0.95 | ? |
| Mut | 1 | muatation speed |
| Imm_dur | 42 | immunity duration |
| Mob | -0.4 (t < 365*1.5) <br> -0.2 (otherwise) | change in mobility |
| Amp | 0.35 | 0.35 seasonal <br> 0.1 endemic |
| pr | 1 - d - l | vaccine protection (<1) |
| pr_ni | pr * 1.2 | recovery protection |
| cont | $Amp*sin(1/180*(pi*2*(t+45)))+1+Mob$ | parameters to change transmission rate |
| naturally infected | 
| beta_naive | beta * cont | $S_0 \Rightarrow I_n$ |
| beta_eff_ni | beta * cont * (1 - pr_ni) | $S_R \Rightarrow I_n$ |
| beta_eff_ni_loss | beta * cont * (1-pr*(1-Loss)) | $S_L \Rightarrow I_n$ |
| vaccinated infected |
| beta_eff | beta*cont*(1-pr) | $V_F \Rightarrow I$ |
| beta_eff_loss | beta*cont*(1-pr*(1-Loss)) | $V_1, V_L \Rightarrow I$ |
| i_total | i + i_n | all infected population |

## 4.	Functions
```
s_dot = -beta_eff_loss*s*i_total - alpha*s + wan*sr; 
```
 $V_L$ = - $V_L$ infected ( $V_L$ $\rightarrow$ $I$ )- $V_L$ take vaccine ( $V_L\rightarrowV_F$ ) + $V_F$ lose immunity ($V_F\rightarrowV_L$)

```
i_dot = beta_eff_loss*s*i_total + beta_eff_loss*v1*i_total + beta_eff*sr*i_total - gamma*i;
```
 $I$ = $V_L, V_1, V_F$ infected ($V_L,V_1,V_F \rightarrow I$)– $I$ recovered ($I \rightarrow R$)

```
sr_dot = alpha*(s + sr + sl_n + sr_n) - beta_eff*sr*i_total - alpha*sr - wan*sr;
```
 $V_F$ = $V_L,V_F,S_L,S_R$ take vaccines – $V_F$ infected - $V_F$ repeated take vaccines - $V_F$ lose immunity ($V_F \rightarrow V_L$)

```
v1_dot = alpha*s0_n - wan*v1 - beta_eff_loss*v1*i_total;
```
 $V_1$ = $S_0$ take first vaccine ($S_0 \rightarrow V_1$) – v1 lose immunity ($V_1 \rightarrow S_0$) – v1 infected

```
s0_n_dot = -beta_naive*s0_n*i_total - alpha*s0_n + wan*v1;
```
$S_0$ = -$S_0$ infected ($S_0 \rightarrow I_n$) – $S_0$ takes first vaccine ($S_0 \rightarrow V_1$) + $V_1$ lose immunity ($V_1 \rightarrow S_0$)

```
i_n_dot = beta_naive*s0_n*i_total + beta_eff_loss_ni*sl_n*i_total + beta_eff_ni*sr_n*i_total - gamma*i_n;
```
$I_n$ = $S_0, S_L, S_R$ infected ($S_0,S_L, S_R \rightarrow I_n$) - $I_n$ recovered ($I_n \rightarrow R$)

```
r_n_dot = gamma*i + gamma*i_n - 1/Imm_dur*r_n;
```
$R$ = $I, I_n$ recovered ($I,I_n \rightarrow R$) - $R$ loses part immunity ($R \rightarrow S_R$)

```
sr_n_dot = 1/Imm_dur*r_n - beta_eff_ni*sr_n*i_total - alpha*sr_n - wan*sr_n;
```
$R$ lose part of immunity ($R \rightarrow S_R$) – $S_R$ infected ($S_R \rightarrow I_n$) – $S_R$ takes vaccines – $S_R$ loses part of immunity ($S_R \rightarrow S_L$)

```
sl_n_dot = -beta_eff_loss_ni*sl_n*i_total - alpha*sl_n + wan*sr_n;
```
$S_L$ = - $S_L$ infected ($S_L \rightarrow I_n$)– $S_L$ takes vaccines + $S_R$ loses part of immunity ($S_R \rightarrow S_L$)


```
a_dot = d*beta_eff_loss*s*i_total + d*beta_eff_loss*v1*i_total + d*beta_eff*sr*i_total + d*beta_naive*s0_n*i_total + d*beta_eff_loss_ni*sl_n*i_total + d*beta_eff_ni*sr_n*i_total + d*alpha*(s+sr+sl_n+sr_n); 
```
$A$ = anigenitic distance * ($V_1, V_F, V_L, S_0, S_L, S_R$ infected + $V_1, V_F, V_L, S_0, S_L, S_R$ takes vaccines)

```
v_dot = Mut*((pr*(1-Loss)*beta_eff_loss*s*i_total + pr*(1-Loss)*beta_eff_loss*v1*i_total + pr*beta_eff*sr*i_total) + 0*beta_naive*s0_n*i_total + (pr_ni*(1-Loss)*beta_eff_loss_ni*sl_n*i_total + pr_ni*beta_eff_ni*sr_n*i_total));
```
$V$ = constant of mutation * ($V_1, V_L, V_F, S_0, S_L, S_R$ different protections)

```
l_dot = 0;
```
Deleterious effect default to be zero

```
ci_dot = beta_eff_loss*v1*i_total + beta_eff_loss*s*i_total + beta_naive*s0_n*i_total + beta_eff_loss_ni*sl_n*i_total + beta_eff*sr*i_total + beta_eff_ni*sr_n*i_total;
```
$CI$ = all infected population ($S_0, V_1, V_L, V_F, S_L, S_R$)

## Figures
Figure.1: natural infected & infected with antibody
<br>
Data: $I$ = x(2), $I_n$ = x(9)

![alt text](./Figure%201.png)

Figure.2: vaccine protected & natural infection protected & full protection
<br>
Data: $S_F$ = x(8),  $S_L$ = x(12), $V_L$ = x(1)

![alt text](./Figure%202.png)

Figure.3: susceptible & natural infected after waning & vaccinated after waning
<br>
Data: $S_O$ = x(8),  $S_L$ = x(12), $V_F$ = x(1)

![alt text](./Figure%203.png)

Figure.4: population antigenicity & virus antigenicity
<br>
Data: $A$ = x(5),  $V$ = x(6)

![alt text](./Figure%204.png)
