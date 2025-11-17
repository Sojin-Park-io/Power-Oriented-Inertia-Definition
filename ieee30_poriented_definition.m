%% === STEP 1: Basic Setup === 
clc; clear; close all
define_constants;
mpc = loadcase('case30');  % Load IEEE 30-bus test case

%% === STEP 2: Data Extraction ===
gen = mpc.gen;
gencost = mpc.gencost;

% Generator power limits
Pg_min = gen(:, PMAX)*0.3;
Pg_max = gen(:, PMAX);

% Load and renewable parameters
load_total = sum(mpc.bus(:, PD));
re_rate = 0;  % Renewable energy rate

% IBR parameters
Pibr_max = 100;  % Max power [MW]
Pibr_min = 0;    % Min power [MW]
eta = 0.9;       % Roundtrip Efficiency
hour = 1;        % Duration time [h]
Eibr_cap = Pibr_max*hour;      % Energy capacity [MWh]
Eibr_min = Eibr_cap*0.2;       % Min state of charge=0.2
Eibr_max = Eibr_cap*0.8;       % Max state of charge=0.8
Eibr_0 = Eibr_cap*0.5;         % Initial state of charge=0.5

% ESS inertia parameters
nibr = 1;
Ribr = 20;        % Droop constant
Hibr_min = 0;     % Min virtual inertia [MW/Hz/s]
Hibr_max = 50;    % Max virtual inertia [MW/Hz/s]

ibr_cost = 3;     % ESS cost [$/MW]

% Generator cost coefficients (quadratic cost function)
a = gencost(:, 5);  % Quadratic coefficient
b = gencost(:, 6);  % Linear coefficient
c = gencost(:, 7);  % Constant term
nGen = size(gen, 1);

%% Initialize result storage
prices_box = [];
H_box = [];
settlement_box = [];
dispatch_box_gen = [];
dispatch_box_ess = [];
cost_box = [];
largest_nonzero_box = {};
largest_box = [];
Energy_box = [];
largest_contingency_box =[];
duals_box = [];
AS_box = [];
duals_P_GEN_box = [];
duals_P_ESS_box = [];
duals_E_ESS_box =[];
nonne_box = [];
duals_sys_box =[];

%% Run optimization for two cases
for j=1:2
    % System inertia and frequency parameters
    H_base = 1;  % Base inertia constant [s] - adjust this to set SI level (e.g., 1.0s for low inertia)
    H = H_base*ones(1,nGen);
    Tpfr = 6;    % PFR required activation time [s]
    T = (1/12);  % Time interval: 5min/60h
    
    alpha = 1;
    beta = 0;
    f0 = 60;     % Nominal frequency [Hz]
    R = 20;      % Droop constant (1/R)
    
    % Frequency stability limits
    RoCoF_lim = 1;      % Rate of Change of Frequency limit [Hz/s]
    Nadir_lim = 0.8;    % Frequency nadir limit [Hz]
    QSS_lim = 0.5;      % Quasi-steady state limit [Hz]
    
    % PFR capabilties
    PFRr = R*Pg_max*Nadir_lim/f0;           % For generators
    PFRd = R*Pg_max*QSS_lim/f0;
    PFRr_ibr = Ribr*Pibr_max*Nadir_lim/f0;  % For ESS
    PFRd_ibr = Ribr*Pibr_max*QSS_lim/f0;

    % Demand setting
    re_rate = 0;
    net_load = load_total*(1-re_rate);

    %% === STEP 3: Build YALMIP Optimization Model ===
    yalmip('clear');
    
    % Decision variables
    Pg = sdpvar(nGen, 1);    % Generator power output
    Pess = sdpvar(nibr, 1);  % ESS power output
    Psi = sdpvar(nGen, 1);   % Synchronous inertia from generators
    Pvi = sdpvar(nibr, 1);   % Virtual inertia from IBR
    Ppfrr = sdpvar(nGen, 1); % PFR reserve (nadir)
    Ppfrr_ibr = sdpvar(nibr, 1);
    Ppfrd = sdpvar(nGen, 1); % PFR reserve (steady-state)
    Ppfrd_ibr = sdpvar(nibr, 1);
    largest_C = sdpvar(1,1); % Largest contingency
    H_ibr = sdpvar(1,1);     % IBR virtual inertia constant
    
    % ESS energy variables
    E_ibr = sdpvar(1,1);     % Final energy state
    Evi = sdpvar(1,1);       % Energy for virtual inertia
    Epfr = sdpvar(1,1);      % Energy for PFR
    
    loss = sdpvar(1,1);
    lossD = sdpvar(1,1);     % Discharge loss
    lossC= sdpvar(1,1);      % Charge loss
    
    % System-level variables
    Intot = sdpvar(1,1);     % Total inertia
    PFRrtot = sdpvar(1,1);   % Total PFR reserve (nadir)
    PFRdtot = sdpvar(1,1);   % Total PFR reserve (steady-state)
    
    % Objective function: minimize total cost
    cost_g = a.* Pg.^2 + b.* Pg + c;
    cost_gen = sum(cost_g); 
    cost_ibr = sum(ibr_cost*(abs(Pess+loss)));
    cost = cost_gen + cost_ibr;
    
    %% Generator constraints
    constraints = [];
 
    Genminmax = Pg_min <= Pg <= Pg_max;              % Power limits
    GenASminmax1 = Pg + Ppfrr <= Pg_max;             % Headroom for PFRs
    GenASminmax2 = Pg + Ppfrd <= Pg_max;
    
    GenIn = Psi == 2*H'.*Pg_max*RoCoF_lim/f0;        % Synchronous inertia
    GenPFRr = Ppfrr <= PFRr;                         % PFR limits
    GenPFRd = Ppfrd <= PFRd;
    GenPFRrd = Ppfrd <= Ppfrr;                       
    
    constraints = [constraints, Genminmax, GenASminmax1, GenASminmax2, ...
                   GenIn, GenPFRd, GenPFRr, GenPFRrd];
    
    %% ESS constraints
    ESSminmax = Pibr_min <= Pess <= Pibr_max;         % Power limits
    ESSASminmax1 = (Pess + alpha*(Pvi+Ppfrr_ibr)) <= Pibr_max;      % Headroom for inertia and PFRs
    ESSASminmax2 = Pibr_min <= (Pess + beta*(Ppfrr_ibr) - Pvi/eta); % Footroom for inertia
    
    ESSIn = Pvi == H_ibr*2*Pibr_max*RoCoF_lim/f0;    % Virtual inertia
    ESSInvar = Hibr_min <= H_ibr <= Hibr_max;
    
    ESSPFRr = Ppfrr_ibr <= PFRr_ibr;                 % PFR limits
    ESSPFRd = Ppfrd_ibr <= PFRd_ibr;
    ESSPFRrd = Ppfrd_ibr <= Ppfrr_ibr;
    
    % Energy constraints
    ESSeminmax = Eibr_min <= E_ibr <= Eibr_max;
    ESSeminmax1 = (Eibr_min + Evi/sqrt(eta) + Epfr/sqrt(eta)) <= Eibr_0;
    ESSeminmax2 = (Eibr_min + Evi/sqrt(eta) + Epfr/sqrt(eta)) <= E_ibr;
    ESSeminmax3 = E_ibr == Eibr_0 - (loss+Pess)*T;  % Energy balance
    
    % Loss calculation
    ESSlossC = lossC == (sqrt(eta)-1)*Pess;
    ESSlossD = lossD == (1/sqrt(eta)-1)*Pess;
    ESSloss1 = lossD <= loss;
    ESSloss2 = lossC <= loss;
    
    % Energy allocation for services
    ESSevi = Evi == Pvi*Nadir_lim/RoCoF_lim/3600;
    ESSepfr = Epfr == (1/2)*Ppfrr_ibr*Tpfr*(1/3600) + Ppfrd_ibr*(T - Tpfr*(1/3600));
    
    % Case 1: positive only, Case 2: bidirectional
    if j == 2
        constraints = [constraints, ESSminmax, ESSASminmax1, ESSASminmax2, ...
            ESSIn, ESSInvar, ESSPFRr, ESSPFRd, ESSPFRrd, ESSeminmax, ...
            ESSeminmax1, ESSeminmax2, ESSeminmax3, ESSlossC, ESSlossD, ...
            ESSloss1, ESSloss2, ESSevi, ESSepfr];
    else
        constraints = [constraints, ESSminmax, ESSASminmax1, ...
            ESSIn, ESSInvar, ESSPFRr, ESSPFRd, ESSPFRrd, ESSeminmax, ...
            ESSeminmax1, ESSeminmax2, ESSeminmax3, ESSlossC, ESSlossD, ...
            ESSloss1, ESSloss2, ESSevi, ESSepfr];
    end
    
    %% System-level constraints
    
    supply_balance = net_load == sum(Pg) + sum(Pess); % Supply-demand balance
    Largen = [Pg <= largest_C; Pess <= largest_C];  % Contingency constraints
    constraints = [constraints, supply_balance, Largen];
    
    % Frequency constraints (RoCOF, QSS)
    RoCoF = largest_C <= Intot;
    QSS = largest_C <= PFRdtot;
    
    % Nadir constraint (second-order cone)
    x = [Intot; PFRrtot; largest_C];
    A = [1/(2*RoCoF_lim), -1/Tpfr, 0;
         0,                 0,        1/sqrt(Nadir_lim)];
    d = [1/(2*RoCoF_lim); 1/Tpfr; 0];
    Nadir = cone([d'*x; A*x]);

    % Aggregate ancillary services
    IN = Intot == sum(Psi)+sum(Pvi);
    PFRR = PFRrtot == sum(Ppfrr)+sum(Ppfrr_ibr);
    PFRD = PFRdtot == sum(Ppfrd)+sum(Ppfrd_ibr);
    
    constraints = [constraints, RoCoF, Nadir, QSS, IN, PFRR, PFRD];
    
    %% Non-negativity constraints
    PGENVAR = [Pg, Psi, Ppfrr, Ppfrd] >= 0;
    PSYSVAR = [largest_C, Intot, PFRrtot, PFRdtot, Pvi, Ppfrr_ibr, ...
               Ppfrd_ibr, H_ibr, Evi, Epfr, E_ibr] >= 0;
    
    constraints = [constraints, PGENVAR, PSYSVAR];
    
    %% === STEP 4: Solve Optimization ===
    options = sdpsettings('solver', 'gurobi', 'gurobi.qcpdual',1);
    diagnostics = optimize(constraints, cost, options);
    error = 10^-4;
    
    if diagnostics.problem ~= 0
        % Infeasible case
        Pg_val = zeros(1, nGen + nibr);
        Pin_val = zeros(1, nGen + nibr);
        Ppfrr_val = zeros(1, nGen + nibr);
        Ppfrd_val = zeros(1, nGen + nibr);
        Largest_val = 0;
        duals = zeros(1,5);
        prices = zeros(1, nGen + nibr + 4);
        settlement = zeros(1, nGen + nibr);
        E_stat_val = zeros(1,5);
    else
        % Extract optimal values
        Pg_val = [value(Pg'), value(Pess)];
        Pin_val = [value(Psi'), value(Pvi)];
        Ppfrr_val = [value(Ppfrr'), value(Ppfrr_ibr)]; 
        Ppfrd_val = [value(Ppfrd'), value(Ppfrd_ibr)];
        
        E_val = [value(Evi), value(Epfr), value(lossC), value(lossD), value(loss)];
        E_stat_val = [value(Eibr_0), value(Evi/sqrt(eta)), value(Epfr/sqrt(eta)), ...
                      value(-(loss+Pess)*T), value(E_ibr)]; 
        
        Largest_val = value(largest_C);

        % Extract dual variables (shadow prices)
        dual_val = dual(supply_balance);
        dual_RoCoF = dual(RoCoF);
        dual_nadir = dual(Nadir);
        dual_QSS = dual(QSS);
        duals = [dual_RoCoF; -dual_nadir(2); -dual_nadir(3); dual_nadir(1); dual_QSS];
        dual_lar = dual(Largen);

        % Individual constraint duals
        dual_indi = [dual(Genminmax); dual(GenASminmax1); dual(GenASminmax2)];
        dual_indi2 = [dual(ESSASminmax2)];
        dual_indi3 = [dual(ESSASminmax1)];
        
        nonne = [dual(PGENVAR); dual(PSYSVAR)];
        duals_sys = [dual(supply_balance); dual(IN); dual(PFRR); dual(PFRD); ...
                     dual_lar; duals];
        
        duals_P = [dual(Genminmax); dual(GenASminmax1); dual(GenASminmax2); ...
                   dual(GenPFRr); dual(GenPFRd); dual(GenPFRrd); dual(GenIn)];
        duals_P_ESS = [dual(ESSminmax); dual(ESSASminmax1); dual(ESSASminmax2); ...
                       dual(ESSPFRr); dual(ESSPFRd); dual(ESSPFRrd); dual(ESSIn); ...
                       dual(ESSInvar)];
        duals_E_ESS = [dual(ESSeminmax1); dual(ESSeminmax2); dual(ESSepfr); ...
                       dual(ESSevi); dual(ESSeminmax3); dual(ESSeminmax); ...
                       dual(ESSloss1); dual(ESSloss2); dual(ESSlossC); dual(ESSlossD)];

        % Calculate market prices for each service
        price_energy = dual_val*ones([nGen+nibr,1]);
        largest_nonzero = find(dual_lar > error);
        largest_nonzero_box{end+1} = largest_nonzero;
        price_energy(largest_nonzero) = dual_val - dual_lar(largest_nonzero);

        price_inertia = dual_RoCoF + (1/2/RoCoF_lim)*(dual_nadir(1)+dual_nadir(2)); %note: signs matter
        price_PFRr = (1/Tpfr)*(-dual_nadir(2)+dual_nadir(1));                       %note: signs matter
        price_PFRd = dual_QSS;
        price_largest = dual_RoCoF + dual_QSS - (1/sqrt(Nadir_lim))*dual_nadir(3);  %note: signs matter
        
        prices = [dual_val, price_inertia, price_PFRr, price_PFRd];
    
        % Calculate settlements (revenue for each resource)
        settlement_G = price_energy(1:6)'.*Pg_val(1:6) + price_inertia*Pin_val(1:6) + ...
                       price_PFRr*Ppfrr_val(1:6) + price_PFRd*Ppfrd_val(1:6);
        settlement_E = price_energy(end)'.*Pg_val(end) + price_inertia*Pin_val(end) + ...
                       price_PFRr*Ppfrr_val(end) + price_PFRd*Ppfrd_val(end);
        settlement = [settlement_G'; settlement_E];
    end
    
    % Store results
    dispatch_box_gen = [dispatch_box_gen, [Pg_val(1:6)'; Pin_val(1:6)'; ...
                        Ppfrr_val(1:6)'; Ppfrd_val(1:6)']];
    dispatch_box_ess = [dispatch_box_ess, [Pg_val(end)'; Pin_val(end)'; ...
                        Ppfrr_val(end)'; Ppfrd_val(end)']];
    Energy_box = [Energy_box, E_stat_val'];

    cost_box = [cost_box, (diagnostics.problem == 0)*value(cost)];
    prices_box = [prices_box, prices'];
    settlement_box = [settlement_box, settlement];

    largest_contingency_box = [largest_contingency_box, Largest_val];
    largest_box = [largest_box, dual_lar];
    AS_box = [AS_box, [sum(Pin_val); sum(Ppfrr_val); sum(Pin_val)/2; sum(Ppfrr_val)/6]];

    duals_box = [duals_box, duals];
    duals_P_GEN_box = [duals_P_GEN_box, duals_P];
    duals_P_ESS_box = [duals_P_ESS_box, duals_P_ESS];
    duals_E_ESS_box = [duals_E_ESS_box, duals_E_ESS];
    nonne_box = [nonne_box, nonne];
    duals_sys_box = [duals_sys_box, duals_sys];
end

% Final results comparison
final_results = [cost_box', prices_box'; 
                 cost_box(1) - cost_box(2), (prices_box(:,1)-prices_box(:,2))'];

%% === Visualization: Power & Energy Allocation ===
nCases = 2;
labels_gen = {'G1','G2','G3','G4','G5','G6','IBR'};
bar_labels_power = {'Inertia^{↓}', 'Energy', 'PFRd', 'PFRr-PFRd', 'Unused', 'Inertia^{↑}'};
bar_labels_energy = {'Energy','Inertia', 'PFR','Unused','Empty'};
incase_spec = ["Positive","Positive","Bidirectional","Bidirectional"];
lccase_spec = ["LargestFixed","LargestVariable","LargestVariable","LargestVariable"];
cost_spec = [0.0, 0.0, 0.0, 13.2];

% Color schemes
colors_power = [
    0.2, 0.6, 1.0;  % Blue
    1.0, 0.4, 0.4;  % Red
    1.0, 0.8, 0.0;  % Yellow
    0.2, 0.8, 0.2;  % Green
    0.6, 0.6, 0.6;  % Gray
    0.2, 0.6, 1.0;  % Blue
];

colors_energy = [
    1.0, 0.4, 0.4;  % Red
    0.2, 0.6, 1.0;  % Blue
    1.0, 0.8, 0.0;  % Yellow
    0.6, 0.6, 0.6;  % Gray
    0.8, 0.8, 0.8;  % Light gray
];

% Create tiled layout
figure('Name', 'Power & Energy Capacity Allocation', ...
       'Position', [100, 100, 300*nCases, 800]);

t = tiledlayout(4, nCases, 'TileSpacing', 'compact', 'Padding', 'compact');

% Store bar handles for legend
h_bar_power = gobjects(1, numel(bar_labels_power));
h_bar_energy = gobjects(1, numel(bar_labels_energy));

for j = 1:nCases
    % Extract power dispatch results
    Pg_final    = dispatch_box_gen(1:nGen, j);
    Ppfrr_final = dispatch_box_gen(nGen*2+1 : nGen*3, j);
    Ppfrd_final = dispatch_box_gen(nGen*3+1 : nGen*4, j);
    Psi_final   = dispatch_box_gen(nGen+1 : nGen*2, j);
    Pess_final  = dispatch_box_ess(1:nibr, j);
    Ppfrr_ibr_final = dispatch_box_ess(nibr*2+1 : nibr*3, j);
    Ppfrd_ibr_final = dispatch_box_ess(nibr*3+1 : nibr*4, j);
    Pvi_final = dispatch_box_ess(nibr+1 : nibr*2, j);

    % Calculate remaining capacity
    remaining_cap = Pg_max - (Pg_final + Ppfrr_final);
    pfr_remaining = Ppfrr_final - Ppfrd_final;

    if Pess_final <= 0
        remaining_cap_ibr = Pibr_max - (Ppfrr_ibr_final + Pvi_final);
        pfr_remaining_ibr = Ppfrr_ibr_final - Ppfrd_ibr_final;
    else
        remaining_cap_ibr = Pibr_max - (Pess_final + Ppfrr_ibr_final + Pvi_final);
        pfr_remaining_ibr = Ppfrr_ibr_final - Ppfrd_ibr_final;
    end

    % Stack power allocation data
    data_stack = [
        -Psi_final, Pg_final, Ppfrd_final, pfr_remaining, remaining_cap, Psi_final;
        0, Pess_final, Ppfrd_ibr_final, pfr_remaining_ibr, remaining_cap_ibr, Pvi_final
    ];

    % Plot power allocation (horizontal bar chart)
    ax1 = nexttile(t, j, [3 1]);
    ax1.YAxis.FontSize = 12;
    ax1.YAxis.FontWeight = 'bold';
    hold(ax1, 'on');
    h = barh(ax1, 1:nGen+nibr, data_stack, 0.6, 'stacked', 'BaseValue', 0);
    for k = 1:length(h)
        h(k).FaceColor = colors_power(k, :);
        h_bar_power(k) = h(k);
    end
    
    % Add capacity limit lines
    for g = 1:nGen
        plot([Pg_max(g), Pg_max(g)], [g-0.4, g+0.5], 'r--', 'LineWidth', 1.75);  
        plot([Pg_min(g), Pg_min(g)], [g-0.4, g+0.5], 'k--', 'LineWidth', 1.75);  
    end
    
    g = nGen + 1;  
    plot([Pibr_max, Pibr_max], [g-0.45, g+0.5], 'r--', 'LineWidth', 1.75);
    plot([Pibr_min, Pibr_min], [g-0.45, g+0.5], 'k--', 'LineWidth', 1.75);

    % Formatting
    xlim([-20 120]);
    xticks(-20:20:120);
    yticks(1:nGen+nibr);
    yticklabels(labels_gen);
    xlabel('Power [MW]','FontSize', 12, 'FontWeight', 'bold');
    if j == 1
        title('(a) Only Positive','FontSize', 12,'FontWeight', 'bold');
    elseif j == 2
        title('(b) Bidirectional','FontSize', 12,'FontWeight', 'bold');
    elseif j == 3
        title ('(c) Increased VI cost','FontSize', 11,'FontWeight', 'bold');
    else
        title ('(d) Increased VI cost & Fixed Contingency','FontSize', 11,'FontWeight', 'bold');
    end
    grid on;

    % Add text annotations
    text(Pess_final/2-10, nGen + 1, sprintf('%.2f', Pess_final), ...
        'FontSize', 10, 'FontWeight', 'bold','Color', 'k', 'VerticalAlignment', 'middle');
    text(sum(data_stack(end,:))-30, nGen + 1, sprintf('%.2f', Pvi_final), ...
        'FontSize', 10, 'FontWeight', 'bold','Color', 'k', 'VerticalAlignment', 'middle');
    text(-4, nGen + 2, sprintf('Largest Contingency = %.2f', largest_contingency_box(:,j)), ...
        'FontSize', 10, 'FontWeight', 'bold','Color', 'k', 'VerticalAlignment', 'middle');
end

% Add legend for power plots
legend(ax1, [h_bar_power, plot(nan, nan, 'r--', 'LineWidth', 2), ...
       plot(nan, nan, 'k--', 'LineWidth', 2)], ...
       [bar_labels_power,'Max Capacity', 'Min Capacity'], ...
       'Location', 'eastoutside', 'Orientation', 'vertical', ...
       'FontSize', 9, 'Box', 'on');

% Plot ESS energy allocation
for j = 1:nCases
    Eibr_0 = Energy_box(1, j);
    Evi_eff = Energy_box(2, j);
    Epfr_eff = Energy_box(3, j);
    Echg = Energy_box(4, j);
    E_final = Energy_box(5, j);

    % Calculate remaining energy capacity
    remaining_up = min(Eibr_cap, Eibr_cap - Echg);
    remaining_dw = min(Eibr_cap, Eibr_cap - Epfr_eff - Evi_eff + Echg);
    remaining_dw = min(remaining_dw, Eibr_cap - Epfr_eff - Evi_eff);

    data_stack_energy = [Echg; -Evi_eff; -Epfr_eff; -remaining_dw; remaining_up];

    % Plot energy allocation
    ax2 = nexttile(t, j + 3*nCases, [1 1]);
    ax2.YAxis.FontSize = 12;
    ax2.YAxis.FontWeight = 'bold';
    hold(ax2, 'on');
    h_e = barh(ax2, 1, data_stack_energy', 'stacked', 'BarWidth', 0.6, 'BaseValue', 0);
    for k = 1:length(h_e)
        h_e(k).FaceColor = colors_energy(k, :);
        h_bar_energy(k) = h_e(k);
    end
    
    % Add SoC limit lines
    h_max = plot([30, 30], [0.6, +1.4], 'r--', 'LineWidth', 1.75);
    h_min = plot([-30, -30], [0.6, +1.4], 'k--', 'LineWidth', 1.75); 
    
    text(30-11, 1.7, 'Max SoC', 'Color', 'r', 'FontSize', 8, 'FontWeight', 'bold');
    text(-30-10, 1.7, 'Min SoC', 'Color', 'k', 'FontSize', 8, 'FontWeight', 'bold');

    xlim([-Eibr_cap/2, Eibr_cap/2]);
    xticks(-50:20:50);
    xticklabels(0:20:100);
    xlabel('SoC [MWh]','FontSize', 12, 'FontWeight', 'bold');
    yticks(1);
    yticklabels('IBR');
    grid on;
end

% Add legend for energy plots
legend(ax2, [h_bar_energy, h_max, h_min], ...
       [bar_labels_energy, 'Max Capacity', 'Min Capacity'], ...
       'Location', 'eastoutside', 'Orientation', 'vertical', ...
       'FontSize', 9, 'Box', 'on');