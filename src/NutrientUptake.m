%% Nutrient uptake model
% Accompaniment to Casey and Follows, 2020; PLOS Computational Biology

% Link to full text: 

% Please refer to the article for details of how parameters were derived or
% measured.

% Play around with the parameters to see how the model behaves with your
% favorite transporter! If you're interested in the nuts and bolts of how
% to get VmaxG using FBA, or in measuring A for your transporter, let me
% know and I'll share some pointers on how that can be done as well. I
% didn't include here because of all the dependencies, but it's fairly
% easy.

% John R. Casey - 05/18/2020
% https://jrcasey.github.io

% Fork or clone from:
% https://github.com/jrcasey/NutrientUptake

% Questions? Feedback?
% jrcasey@mit.edu

%% Compile parameters for E. coli growth on glucose

VmaxG = 1.77e-18; % Maximum uptake rate in replete batch (mol cell-1 s-1), derived from FBA
nG = 2775; % Transporter abundance in replete batch (cell-1), measured by Schmidt et al., 2016
r = 8.785e-7; % Cell radius (m), measured by Schmidt et al., 2016
A = 3.9147e-17; % PtsG catchment area (m2), measured using modeled protein structure.
T = 37; % Temperature (C), measured by Schmidt et al., 2016
MW = 180.156; % Molecular weight of glucose (Da). Be sure that your molecule is hydrated before calculating MW
u = 0; % Advective velocity (m s-1)

%% define a vector of substrate concentrations (easier to see in log space)

S_res = 500; % resolution in the S domain
S = logspace(-4,2,S_res); % Nutrient concentration range (mole m-3). You may need to adjust this for another transport process

%% Compute molecular diffusivity

D = getDiffusivity(MW,T); % Hydrated molecular diffusivity (m2 s-1) 

%% Compute transporter in vivo catalytic rate

kcat = VmaxG./nG; % catalytic constant (molecules tp-1 s-1)

%% Compute maximum number of transporters

f_max = 0.085; % Maximum fraction of the cell surface that can be covered by transporters... wonkiest term in here. See text for a discussion of n_max.
n_max = f_max.*(4.*pi().*(r.^2))./A; % cell-1

%% Get glucose-replete batch-acclimated parameters

[ks_d, ks_p] = getK(VmaxG, kcat, nG, r, A, D, u); % porter limit and diffusive limit of half saturation concentration (mol m-3)

%% Calculate various critical concentrations

% Lower bound to growth limited domain of n-star
S_G_lb = (nG.*kcat) ./ (D.*r); % mol m-3
S_G_lb_idx = min(find(S-S_G_lb >= 0)); % nearest concentration index in S

% Calculate S-star. This is the positive root of the polynomial
% intersection of n_optP and n_optD.
kprime = (kcat./(r .* D));
a = 1./(kprime.*nG);
b = -2;
c = -ks_p;
p = [a b c];
y0 = roots(p);
S_star = y0(1); % mol m-3
S_star_idx = min(find(S-S_star >= 0)); % nearest concentration index in S

% Lower and upper bounds of surface area limitation domain (true if optimal
% n excedes n_max
max_n_opt = (ks_p + S(S_star_idx)) ./ ((S(S_star_idx)./nG) - (kcat./(D.*r)));
if max_n_opt > n_max
    S_SA_lb = (n_max .* kcat) ./ (r .* D); % mol m-3
    S_SA_lb_idx = min(find(S-S_SA_lb >= 0)); % nearest concentration index in S
    S_SA_ub = ( ((n_max .* kcat)./(D.*r)) + ks_p ) ./ ( (n_max./nG) -1 ); % mol m-3
    S_SA_ub_idx = min(find(S-S_SA_ub >= 0)); % nearest concentration index in S
    SA_transition = true;
else
    SA_transition = false;
end

%% Compute optimal number of transporters

% In the growth limit range, such that v = vmax
n_optG = (ks_p + S) ./ ((S./nG) - (kcat./(D.*r))); % cell-1
% Constrain to the feasible domain (nans outside)
n_optG(1:S_G_lb_idx-1) = NaN;

% In the diffusive domain, such that v = vD
n_optD = (S.*r.*D)./kcat; % cell-1

% In the surface area limited domain, such that v = n_max*kcat*S/(ks_p + ks_d + S),
% where ks_d is evaluated at n_max
n_optSA = nan(1,numel(S));
if SA_transition
    n_optSA(S_SA_lb_idx:S_SA_ub_idx) = n_max; % cell-1
end

% Minimum of all three sets
n_opt_minDP = nanmin([n_optG;n_optSA;n_optD]); % cell-1

%% Compute uptake for optimal n

for a = 1:numel(S)
    Vmax(a) = n_opt_minDP(a).*kcat; % mole cell-1 s-1
    [ks_d(a), junk] = getK(Vmax(a), kcat, n_opt_minDP(a), r, A, D, u); % mol m-3
    [v(a)] = getUptake(S(a),Vmax(a),ks_d(a), ks_p); % mole cell-1 s-1
end
ks = ks_d + ks_p; % mol m-3
v_P = v.*1e15.*3600; % fmol cell-1 h-1

%% Michaelis-Menten kinetics for glucose-replete acclimated cells, just for comparison

[ks_d3, junk] = getK(VmaxG, kcat, nG, r, A, D, u); % mol m-3
for a = 1:numel(S)
    [v3(a)] = getUptake(S(a),VmaxG,ks_d3, ks_p); % mole cell-1 s-1
end
ks3 = ks_d3 + ks_p; % mol m-3
v_P3 = v3.*1e15.*3600; % mole cell-1 s-1

%% Compute maximum diffusive flux toward the cell
% Based on Berg and Purcell, 1977
vD = 4*pi()*r .* D .* S .* ( (n_optD.*A) ./ (4.*pi().*(r.^2)) ); % mole cell-1 s-1

%% Uptake in the n vs S plane
n_res = 100; % resolution in the n domain
nLim = 0.6*n_max; % choose some maximum value for the domain... you'll need to adjust this to zoom in on the region you're interested in.
n_Plane = linspace(1,nLim,n_res);
% loop through the plane and compute uptake rate at each pixel
for a = 1:numel(n_Plane)
    Vmax_Plane = n_Plane(a).*kcat; % mole cell-1 s-1
    [ks_d_Plane, ks_p_Plane] = getK(Vmax_Plane, kcat, n_Plane(a), r, A, D, u); % mol m-3
    for b = 1:numel(S)
        [v_Plane(a,b)] = getUptake(S(b), Vmax_Plane, ks_d_Plane, ks_p_Plane); % mole cell-1 s-1
    end
end
v_Plane2 = v_Plane .* 1e15 .* 3600; % fmol cell-1 h-1

%% Add experimental data from Schmidt et al., 2016

nE = [2381 3919 6753 9402]; % transporters per cell corresponding to 0.12, 0.20, 0.35, and 0.50 h-1
rE = 1e-7.*[7.683 7.919 8.306 8.628]; % radius (m)
muE = [0.12 0.20 0.35 0.50]; % growth rate (h-1)
mu_max = 1.80; % maximum growth rate (h-1) using Schmidt's LB batch growth rate

% calculate steady-state glucose concentrations in the chemostats
for a = 1:numel(nE)    
    VmaxE(a) = kcat.*nE(a); % mole cell-1 s-1
    [ks_d_E(a), ks_p_E(a)] = getK(VmaxE(a), kcat, nE(a), rE(a), A, D, u); % mol m-3 
end
ks_E = ks_p_E + ks_d_E; % mol m-3
S_E = (muE .* ks_E) ./ (mu_max - muE); % mol m-3
for a = 1:numel(nE)
    [v_E(a)] = getUptake(S_E(a), VmaxE(a), ks_d_E(a), ks_p_E(a)); % mole cell-1 s-1
end
v_E2 = v_E .* 1e15 .* 3600; % fmol cell-1 h-1

%% Plot results

figure
subplot(2,1,1) % plot uptake in Tp vs S space. Overlay data
h1 = contourf(S,n_Plane,v_Plane2,20,'LineStyle','none');
hold on
h2 = plot(S,n_opt_minDP,'-r','LineWidth',4);
if SA_transition
    h5 = plot([S_SA_lb S_SA_lb],[0 nLim],'--g','LineWidth',3)
    h6 = plot([S_SA_ub S_SA_ub],[0 nLim],'--g','LineWidth',3)
else
    h7 = plot([S_star S_star],[0 nLim],'--g','LineWidth',3)
end
h3 = plot([min(S) max(S)],[n_max n_max],'--k','LineWidth',3)
h4 = plot([S_E 33],[nE nG],'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',15,'LineWidth',3);
xlabel('S_\infty [mol m^-^3]');
ylabel('PtsG [cell^-^1]');
xlim([1e-4 1e2])
ylim([0 nLim])
h = colorbar;
colormap('jet')
ylabel(h, 'Uptake rate [fmol cell^-^1 h^-^1]');
set(gca,'FontSize',20);
set(gca,'XScale','log');
if max(n_opt_minDP) > nLim % Ugly AF but works
    if SA_transition
        hl = legend([h4, h2, h3, h5, h6],[{'Schmidt et al., 2016'},{'n^*'},{'n_m_a_x'},{'S_S_A^l^b'},{'S_S_A^u^b'}])
    else
        hl = legend([h4, h2, h3, h7],[{'Schmidt et al., 2016'},{'n^*'},{'n_m_a_x'},{'S^*'}])
    end
else
    if SA_transition
        hl = legend([h4, h2, h5, h6],[{'Schmidt et al., 2016'},{'n^*'},{'S_S_A^l^b'},{'S_S_A^u^b'}])
    else
        hl = legend([h4, h2, h7],[{'Schmidt et al., 2016'},{'n^*'},{'S^*'}])
    end
end
set(hl,'FontSize',20,'EdgeColor','w')

subplot(2,1,2) % Uptake versus concentration
plot([S_E 33],[v_E2 VmaxG.*1e15.*3600],'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',15,'LineWidth',3);
hold on
plot(S,v_P,'-r','LineWidth',4);
plot(S,v_P3,'-c','LineWidth',4);
plot(S,vD.*1e15.*3600,'--b','LineWidth',3);
plot([S_star S_star],[0 2.*VmaxG.*1e15.*3600],'--g','LineWidth',3)
plot([min(S) max(S)],[VmaxG.*1e15.*3600 VmaxG.*1e15.*3600],'--','Color',[253, 187, 132]./256,'LineWidth',3);
set(gca,'XScale','log');
xlim([1e-4 1e2])
ylim([0 1.2.*VmaxG.*1e15.*3600])
xlabel('S_\infty [mol m^-^3]')
ylabel('Uptake Rate [fmol cell^-^1 h^-^1]')
hl = legend('Schmidt et al., 2016','Optimal','Batch-acclimated','Maximum Diffusive Flux','S^*','V_m_a_x^G','Location','SouthEast')
set(hl,'EdgeColor','w')
set(gca,'FontSize',20)



