function [prdData, info] = predict_Sparus_aurata(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

  if E_Hh >= E_Hb | E_Hh < 0
    info = 0; prdData = []; return
  end

  % compute temperature correction factors
  TC_ah = tempcorr(temp.ah, T_ref, T_A);
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC_aj = tempcorr(temp.aj, T_ref, T_A);
  TC_ap = tempcorr(temp.ap, T_ref, T_A);
  TC_am = tempcorr(temp.am, T_ref, T_A);
  TC_GSI = tempcorr(temp.GSI, T_ref, T_A);
  TC_tL_larv = tempcorr(temp.tL_larv, T_ref, T_A);
  TC_Bav = tempcorr(temp.tL_Bav, T_ref, T_A);
  
  % zero-variate data

  % life cycle
  pars_tj = [g k l_T v_Hb v_Hj v_Hp];
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
  
  % initial
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve

  % hatch   
  [U_H aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
  a_h = aUL(2,1); aT_h = a_h/ TC_ah; % d, age at hatch at f and T
  Lw_h = aUL(2,3)/ del_M_larv;            % cm, total length at hatch
  
  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Lw_b = L_b/ del_M_larv;                % cm, physical length at birth at f
  Ww_b = L_b^3 * (1 + f * w);       % g, wet weight at birth at f 
  aT_b = t_b/ k_M/ TC_ab;           % d, age at birth at f and T

  % metam
  aT_j = t_j/ k_M/ TC_aj;           % d, age at metam
  L_j = L_m * l_j;                  % cm, structural length at metam
  Ww_j = L_j^3 * (1 + f * w);       % g, wet weight at metam
  s_M = L_j/ L_b;                   % -, acceleration
  
  % puberty 
  L_p = L_m * l_p;                  % cm, structural length at puberty at f
  Lw_p = L_p/ del_M_ad;                % cm, physical length at puberty at f
  Ww_p = L_p^3 *(1 + f * w);        % g, wet weight at puberty 
  aT_p = t_p/ k_M/ TC_ap;           % d, age at puberty at f and T

  % ultimate
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M_ad;                % cm, ultimate physical length at f
  Ww_i = L_i^3 * (1 + f * w);       % g, ultimate wet weight 
 
  % reproduction
  GSI = TC_GSI * 365 * k_M * g/ f^3/ (f + kap * g * y_V_E);
  GSI = GSI * ((1 - kap) * f^3 - ...
     k_M^2 * g^2 * k_J * U_Hp/ v^2/ s_M^3); % mol E_R/ mol W

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T
  
  % pack to output
  prdData.ah = aT_h;
  prdData.ab = aT_b;
  prdData.aj = aT_j;
  prdData.ap = aT_p;
  prdData.am = aT_m;
  prdData.Lh = Lw_h;
  prdData.Lb = Lw_b;
  prdData.Lp = Lw_p;
  prdData.Li = Lw_i;
  prdData.Wwb = Ww_b;
  prdData.Wwj = Ww_j;
  prdData.Wwp = Ww_p;
  prdData.Wwi = Ww_i;
  prdData.GSI = GSI;
  
  % uni-variate data
  
  % time-length
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_tL_larv);
  kT_M = k_M * TC_tL_larv; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M;        
  L_b = L_m * l_b;  L_j = L_m * l_j; L_i = L_m * l_i; tT_j = (t_j - t_b)/ kT_M;
  L_bj = L_b * exp(tL_larv((tL_larv(:,1) <= tT_j),1) * rT_j/ 3);
  L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tL_larv((tL_larv(:,1) > tT_j),1) - tT_j)); % cm, expected length at time
  L = [L_bj; L_jm];
  ELw_larv = L ./ del_M_larv;

%   w1 = min(1, max(0, (L - L_b)/ (L_j - L_b))); % smooth transition from larvea to juvenile (metamorphosis)
%   ELw_larv = L ./ (del_M * (1 - w1) + del_M * w1); % cm, total length ?? (in the original file, but w1 does not have an impact)
%   ELw_larv = L ./ (del_M * w1); % cm, total length
  
  
  %% For data added at DEB School 2021
  
% t-L & t-Ww for tLW_Bav
  F = f_Bav;
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, F);
  L_i = l_i * L_m;
  rT_B = rho_B * k_M * TC_Bav; % d, 1/von Bert growth rate, because only adults
  ELw_Bav = (L_i - (L_i - L0.tL_Bav * del_M_ad) * exp( - rT_B * tL_Bav(:,1))) ./ del_M_ad; % cm, total length
  EWw_Bav = (ELw_Bav * del_M_ad).^3 * (1 + F * w); % g, wet weight
  
% L-Ww for LW_Bav
  ELWw_Bav = (LWw_Bav(:,1) * del_M_ad).^3 * (1 + F * w); % g, wet weight


  % pack to output
  prdData.tL_larv = ELw_larv;
  prdData.tL_Bav = ELw_Bav;
  prdData.tWw_Bav = EWw_Bav;
  prdData.LWw_Bav = ELWw_Bav;

