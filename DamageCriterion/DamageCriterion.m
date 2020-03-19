function Fun_crit = damageCriterion(crit_Parameter)


beta        = crit_Parameter.beta;
beta_angle  = crit_Parameter.beta_angle;
hr          = crit_Parameter.hr;
D           = crit_Parameter.Damage;
index1      = abs(1-D)>eps;
index2      = abs(1-D)<eps;
gD          = zeros(length(D),1);
gD(index1)  = (1-D(index1)).*(1-beta*log(1-D(index1)));
gD(index2 ) = 0;

hD          = hr+(1-D).^beta_angle*(1-hr);
c           = crit_Parameter.c;
sR          = crit_Parameter.sR;
phi         = crit_Parameter.phi/180*pi;
tau_c       = (c^2+(sR*tan(phi))^2 ) / (2*sR*tan(phi));

tau         = crit_Parameter.tau;
sn          = crit_Parameter.sn;
Fun_crit = (tau.^2 -(hD.*sn*tan(phi)).^2 + 2*hD.*gD.*tau_c.*sn*tan(phi) -gD.^2*c^2)/tau_c^2;

end