function dPop=SEIR_model(t,pop)

fixed_param;

dPop=zeros(4*n,1);

S=pop(1:4:(4*n)); 
E=pop(2:4:(4*n));
I=pop(3:4:(4*n)); 
R=pop(4:4:(4*n));


%% Movement matrix
% the element (i,j) denote rate from patch j to patch i  ...
% Fully connected
m12=epsilon/((d11^eta1)*nu1(2)); m13=epsilon/((d12^eta1)*nu1(3));
m14=epsilon/((d13^eta1)*nu1(4)); m15=epsilon/((d14^eta1)*nu1(5));
m21=epsilon/((d11^eta1)*nu1(1)); m23=epsilon/((d15^eta1)*nu1(3));
m24=epsilon/((d16^eta1)*nu1(4)); m25=epsilon/((d17^eta1)*nu1(5));
m31=epsilon/((d12^eta1)*nu1(1)); m32=epsilon/((d15^eta1)*nu1(2));
m34=epsilon/((d18^eta1)*nu1(4)); m35=epsilon/((d19^eta1)*nu1(5));
m41=epsilon/((d13^eta1)*nu1(1)); m42=epsilon/((d16^eta1)*nu1(2));
m43=epsilon/((d18^eta1)*nu1(3)); m45=epsilon/((d10^eta1)*nu1(5));
m51=epsilon/((d14^eta1)*nu1(1)); m52=epsilon/((d17^eta1)*nu1(2));
m53=epsilon/((d19^eta1)*nu1(3)); m54=epsilon/((d10^eta1)*nu1(4));

%Ring of patches
% m12=0; m13=0;
% m14=0; m15=epsilon/((d14^eta1)*nu1(5));
% m21=epsilon/((d11^eta1)*nu1(1)); m23=0;
% m24=0; m25=0;
% m31=0; m32=epsilon/((d15^eta1)*nu1(2));
% m34=0; m35=0;
% m41=0; m42=0;
% m43=epsilon/((d18^eta1)*nu1(3)); m45=0;
% m51=0; m52=0;
% m53=0; m54=epsilon/((d10^eta1)*nu1(4));

% Star network
% m12=epsilon/((d11^eta1)*nu1(2)); m13=epsilon/((d12^eta1)*nu1(3));
% m14=epsilon/((d13^eta1)*nu1(4)); m15=epsilon/((d14^eta1)*nu1(5));
% m21=epsilon/((d11^eta1)*nu1(1)); m23=0;
% m24=0; m25=0;
% m31=epsilon/((d12^eta1)*nu1(1)); m32=0;
% m34=0; m35=0;
% m41=epsilon/((d13^eta1)*nu1(1)); m42=0;
% m43=0; m45=0;
% m51=epsilon/((d14^eta1)*nu1(1)); m52=0;
% m53=0; m54=0;


Mig=[0 m12 m13 m14 m15 ; 
     m21 0 m23 m24 m25;
     m31 m32 0 m34 m35;
     m41 m42 m43 0 m45;
     m51 m52 m53 m54 0];
 


dPop(1:4:(4*n))= Pi_h -((beta_h.*I)./(N)).*S - mu_h.*S - sum(Mig)'.*S + Mig*S - Mig*(((alpha_h.*I)./(N)).*S) ;

dPop(2:4:(4*n))=((beta_h.*I)./(N)).*S - (mu_h + gamma_h).*E - sum(Mig)'.*E + Mig*(((alpha_h.*I)./(N)).*S) + Mig*((1-xi_h)*E);

dPop(3:4:(4*n))= gamma_h.*E  - (sigma_h + mu_h + delta_h).*I - sum(Mig)'.*I + Mig*(xi_h*E) + Mig*((1-p_h)*I);

dPop(4:4:(4*n))= sigma_h.*I  - mu_h.*R - sum(Mig)'.*R + Mig*R + Mig*(p_h*I);

end