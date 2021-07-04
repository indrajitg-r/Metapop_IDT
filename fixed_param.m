%% Fixed parameters of the SEIR model
n=5; %%% Number of patches

%%Fixed parameters
mu_h = 1/(52*60);
beta_h=[0.3; 0.24; 0.18; 0.12; 0.06];
alpha_h=[0.4; 0.35; 0.3; 0.25; 0.2];
xi_h = 0.3;
p_h = 0.001;
gamma_h = 0.15;
sigma_h = 0.09; 
delta_h = 0.05; 
d11=100; d12=110; d13=120; d14=130; d15=140;
d16=150; d17=160; d18=170; d19=180; d10=190;
eta1=0.25;
epsilon=0.005;
nu1=[7; 6; 5; 4; 3];


S0=[100000;100000;100000;100000;100000];%% total population at risk
E0=[100; 0; 0; 0; 0]; %% #of exposed at week 1
I0=[10; 0; 0; 0; 0]; %% #of infected at week 1
R0 =[0; 0; 0; 0; 0]; %% #of recovered at week 1


y0(1:4:4*n)=S0;
y0(2:4:4*n)=E0;
y0(3:4:4*n)=I0;
y0(4:4:4*n)=R0;

%%%%% Total population of the patches
N1=S0(1)+E0(1)+I0(1)+R0(1);
N2=S0(2)+E0(2)+I0(2)+R0(2);
N3=S0(3)+E0(3)+I0(3)+R0(3);
N4=S0(4)+E0(4)+I0(4)+R0(4);
N5=S0(5)+E0(5)+I0(5)+R0(5);
N=[N1; N2; N3; N4; N5];

Pi_h= mu_h.*N;