function [Vs, Vd, Va, t, I_output] = doublet_3comp_model(I, inj_site, method, param, T, step_length)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A three-compartment motoneuron model (dendrite + soma + axon) with
% holding + step current protocol tuned to produce doublets
% Lumped compartments: dendritic, soma+AIS, axon (first node of Ranvier)
%
% The model (inspired by single-compartment model by Purvis & Butera 2005)
% was developed for: An axonal sodium current gates motoneuron doublets and
% force amplification (https://doi.org/10.64898/2026.07.29.741588)
%
% Implemented by: Robin Rohlén (last updated 2026-08-16)
%
% INPUTS:
%   I           = [I_hold, I_step] in nA (or model units)
%   inj_site    = injection site 'soma' or 'dendrite'
%   method      = 'euler' or 'RK4'
%   param       = parameter vector
%                   param(1)    gNaP_s      (soma persistent Na)
%                   param(2)    gSK_s       (soma SK)
%                   param(3)    EK          (K reversal potential)
%                   param(4)    g_c_sa      (soma-axon coupling)
%                   param(5)    gNaP_ax     (axon persistent Na)
%                   param(6)    gT_s        (soma T-type Ca)
%                   param(7)    gP_s        (soma P-type Ca)
%                   param(8)    gN_s        (soma N-type Ca)
%                   param(9)    g_c_ds      (dendrite-soma coupling)
%                   param(10)   gT_d        (dendrite T-type Ca)
%                   param(11)   gN_d        (dendrite N-type Ca)
%                   param(12)   tau_hNaP_s  (soma NaP inactivation time)
%                   param(13)   tau_hNaP_ax (axon NaP inactivation time)
%   T           = length of simulation in ms (default 2000 ms)
%   step_length = length of step current (default T-1000 ms = 1000 ms)
%
% OUTPUTS:
%   Vs          = Membrane voltage - somatic compartment (mV)
%   Vd          = Membrane voltage - dendritic compartment (mV)
%   Va          = Membrane voltage - axonal compartment (mV)
%   t           = time vector
%   I_output    = struct with outputs, e.g., currents and kinetics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1, error('Must input a current value.'); end
if nargin < 2, inj_site = 'soma'; end
if nargin < 3, method = 'euler'; end
if nargin < 4, param = [0.03 0.30 -90 0.20 0.04 0.05 0.02 0.02 0.07 0.08 0.15 150 25]; end
if nargin < 5, T=2000; end
if nargin < 6, step_length=T-1000; end

% Add 1.5 s holding current to ensure steady state before step current
T = T + 1500;

if 2500+step_length>T
    error('Step length is too long (shorten or extend T)');
end

dt = 0.001;

% Build stimulation waveform - applied to either soma or dendrite
if length(I)==2
    Nt = round(T/dt);
    t  = (0:Nt-1)'*dt;

    I_stim = I(1)*ones(Nt,1);
    I_stim(round(2500/dt):round((2500+step_length)/dt)) = I(2);
else
    I_stim=I;
    Nt = length(I_stim);
    t  = (0:Nt-1)'*dt;
end

% Initialise state vector (29 states)
y = zeros(Nt, 29);

% Initial conditions
y(1,:) = [-73.8136, -73.8136, -73.8136, ... % Vd, Vs, Va
    0.015, 0.981, 0.002, 0.797, 0.158, ...  % m, h, mNaP, hNaP, n
    0.001, 0.562, 0, 0.001, 0.649, ...      % mT, hT, mP, mN, hN
    0, 0.057, 0.287, 0.182, 0.10, ...       % zSK_s, mA, hA, mH, Ca_s
    0.001, 0.562, 0.001, 0.649, ...         % mT_d, hT_d, mN_d, hN_d
    0, 0.10, ...                            % zSK_d, Ca_d
    0.981, 0.002, ...                       % mNaP_ax, hNaP_ax
    0.01, 0.98, 0.02];                      % m_ax, h_ax, n_ax

% Determine which compartment receives the injection
inject_to_soma = strcmpi(inj_site, 'soma');
inject_to_dend = strcmpi(inj_site, 'dendrite') || strcmpi(inj_site, 'dend');

if ~(inject_to_soma || inject_to_dend)
    error('inj_site must be ''soma'' or ''dendrite''');
end

% Main integration loop
for i = 1:Nt-1
    if inject_to_soma
        Iext_s = I_stim(i);
        Iext_d = 0;
    else
        Iext_s = 0;
        Iext_d = I_stim(i);
    end

    if strcmp(method,'RK4')
        k1 = doublet_3comp_model_rhs(y(i,:), Iext_s, Iext_d, param);
        k2 = doublet_3comp_model_rhs(y(i,:) + dt/2*k1', Iext_s, Iext_d, param);
        k3 = doublet_3comp_model_rhs(y(i,:) + dt/2*k2', Iext_s, Iext_d, param);
        k4 = doublet_3comp_model_rhs(y(i,:) + dt*k3',   Iext_s, Iext_d, param);

        y(i+1,:) = y(i,:) + dt/6*(k1' + 2*k2' + 2*k3' + k4');
    elseif strcmp(method,'euler')
        dydt = doublet_3comp_model_rhs(y(i,:), Iext_s, Iext_d, param);
        y(i+1,:) = y(i,:) + dt*dydt';
    else
        error('No such method.')
    end

    y(i+1, 4:17) = min(max(y(i+1, 4:17), 0), 1); % soma gating + zSK_s
    y(i+1, 19:23) = min(max(y(i+1, 19:23), 0), 1); % dendrite gating + zSK_d
    y(i+1, 25:29) = min(max(y(i+1, 25:29), 0), 1); % axon gating
    y(i+1, 18) = max(y(i+1, 18), 0); % Ca_s
    y(i+1, 24) = max(y(i+1, 24), 0); % Ca_d
end

% Discard pre-stimulus interval
discard_idx = 1:round(1500/dt);
y(discard_idx, :) = [];
t = t(1:size(y,1));
I_stim(discard_idx, :) = [];

% Unpack voltages
Vd = y(:,1);
Vs = y(:,2);
Va = y(:,3);

%% Recompute currents for output (using same params as RHS)

% Soma states
m    = y(:,4);
h    = y(:,5);
mNaP = y(:,6);
hNaP = y(:,7);
n    = y(:,8);
mT   = y(:,9);
hT   = y(:,10);
mP   = y(:,11);
mN   = y(:,12);
hN   = y(:,13);
zSK_s = y(:,14);
mA   = y(:,15);
hA   = y(:,16);
mH   = y(:,17);
Ca_s = y(:,18);

% Dendrite states
mT_d   = y(:,19);
hT_d   = y(:,20);
mN_d   = y(:,21);
hN_d   = y(:,22);
zSK_d  = y(:,23);
Ca_d   = y(:,24);

% Axon states
mNaP_ax = y(:,25);
hNaP_ax = y(:,26);
m_ax    = y(:,27);
h_ax    = y(:,28);
n_ax    = y(:,29);

% Parameters (for current recomputation)
gNaP_s  = param(1); 
gSK_s   = param(2);
EK      = param(3);
g_c_sa  = param(4);
gNaP_ax = param(5); 
gT_s    = param(6);
gP_s    = param(7);
gN_s    = param(8);
g_c_ds  = param(9);
gT_d    = param(10);
gN_d    = param(11);
tau_hNaP_s = param(12);
tau_hNaP_ax = param(13);

gNa   = 1.3;
gK    = 1.5;
gleak = 0.008;
gA    = 1.1;
gH    = 0.005;
Cm_s  = 0.08;
gSK_d   = 0.30;
gleak_d = 0.025;
gNa_ax  = 0.20;
gK_ax   = 0.60; 
gleak_ax = 0.004;

ENa = 60;
Eleak = -70;
ECa = 50;
EH = -38.8;

% Soma currents
INa     = gNa  .* m.^3 .* h .* (Vs-ENa);
INaP_s  = gNaP_s .* mNaP .* hNaP .* (Vs-ENa);
IK      = gK   .* n.^4 .* (Vs-EK);
Ileak_s = gleak .* (Vs-Eleak);
IT_s    = gT_s .* mT .* hT .* (Vs-ECa);
IN_s    = gN_s .* mN .* hN .* (Vs-ECa);
IP_s    = gP_s .* mP .* (Vs-ECa);
ICa_s   = IT_s + IN_s + IP_s;
IA      = gA   .* mA .* hA .* (Vs-EK);
IH      = gH   .* mH .* (Vs-EH);
ISK_s   = gSK_s .* zSK_s.^2 .* (Vs-EK);

% Dendrite currents
Ileak_d  = gleak_d .* (Vd-Eleak);
IT_d     = gT_d .* mT_d .* hT_d .* (Vd-ECa);
IN_d     = gN_d .* mN_d .* hN_d .* (Vd-ECa);
ICa_d    = IT_d + IN_d;
ISK_d    = gSK_d .* zSK_d.^2 .* (Vd-EK);

% Axon currents
INa_ax   = gNa_ax .* m_ax.^3 .* h_ax .* (Va-ENa);
IK_ax    = gK_ax  .* n_ax.^4 .* (Va-EK);
INaP_ax  = gNaP_ax .* mNaP_ax .* (Va-ENa);
Ileak_ax = gleak_ax .* (Va-Eleak);

% Coupling currents (signed for inspection)
I_c_s_to_a = g_c_sa .* (Vs - Va);
I_c_s_to_d = g_c_ds .* (Vs - Vd);

%% Store outputs

% Soma
I_output.INa = INa;
I_output.INaP_s = INaP_s;
I_output.mNaP = mNaP;
I_output.hNaP = hNaP;
I_output.IK   = IK;
I_output.IA   = IA;
I_output.ISK_s = ISK_s;
I_output.IH   = IH;
I_output.Ileak_s = Ileak_s;
I_output.ICa_s = ICa_s;
I_output.IT_s = IT_s;
I_output.IN_s = IN_s;
I_output.IP_s = IP_s;
I_output.Ca_s = Ca_s;

% Dendrite
I_output.Ileak_d = Ileak_d;
I_output.ICa_d = ICa_d;
I_output.IT_d = IT_d;
I_output.IN_d = IN_d;
I_output.ISK_d = ISK_d;
I_output.Ca_d = Ca_d;
I_output.zSK_d = zSK_d;

% Axon
I_output.INa_ax  = INa_ax;
I_output.IK_ax   = IK_ax;
I_output.INaP_ax = INaP_ax;
I_output.mNaP_ax = mNaP_ax;
I_output.hNaP_ax = hNaP_ax;
I_output.Ileak_ax = Ileak_ax;

% Coupling and stim
I_output.I_c_s_to_a = I_c_s_to_a;
I_output.I_c_s_to_d = I_c_s_to_d;
I_output.I_stim = I_stim;
I_output.inj_site = inj_site;

end

function dydt = doublet_3comp_model_rhs(y, Iext_s, Iext_d, param)
% Three-compartment RHS (dendrite + soma + axon)

%% Unpack voltages
Vd = y(1);
Vs = y(2);
Va = y(3);

%% Soma states
m = y(4);
h = y(5);
mNaP = y(6);
hNaP = y(7);
n = y(8);
mT = y(9); 
hT = y(10);
mP = y(11);
mN = y(12); 
hN = y(13);
zSK_s = y(14);
mA = y(15); 
hA = y(16);
mH = y(17);
Ca_s = y(18);

%% Dendrite states
mT_d = y(19);   
hT_d = y(20);
mN_d = y(21);   
hN_d = y(22);
zSK_d = y(23);
Ca_d = y(24);

%% Axon states
mNaP_ax = y(25); 
hNaP_ax = y(26);
m_ax = y(27);   
h_ax = y(28);
n_ax = y(29);

%% Parameters
gNaP_s  = param(1);
gSK_s   = param(2);
EK      = param(3);
g_c_sa  = param(4);
gNaP_ax = param(5);
gT_s    = param(6);
gP_s    = param(7);
gN_s    = param(8);
g_c_ds  = param(9);
gT_d    = param(10);
gN_d    = param(11);
tau_hNaP_s = param(12);
tau_hNaP_ax = param(13);

%% Fixed soma conductances (from pab_rhs_RK4_2comp_Na_K.m)
gNa   = 1.3;
gK    = 1.5;
gA    = 1.1;
gH    = 0.005;
Cm_s  = 0.08;
gleak = 0.008;

gNa_ax  = 0.20;
gK_ax   = 0.60;
Cm_a    = 0.040;
gleak_ax= 0.004;

gSK_d   = 0.30;
Cm_d    = 0.20;
gleak_d = 0.025;

%% If all params are zero, make the cell passive (set all conduct. to 0)
if nnz(param)==0
    gNaP_s  = 0;
    gSK_s   = 0;
    EK      = -90;
    g_c_sa  = 0.20;
    gNaP_ax = 0;
    gNa_ax  = 0;
    gK_ax   = 0;
    gT_s    = 0;
    gP_s    = 0;
    gN_s    = 0;
    gSK_d   = 0;
    g_c_ds  = 0.07;
    gT_d    = 0;
    gN_d    = 0;
    gNa     = 0;
    gK      = 0;
    gA      = 0;
    gH      = 0;
end

%% Ca handling (split pools, same constants as before)
K1 = -0.00053;
K2 = 0.015;

%% Reversal potentials
ENa = 60;
Eleak = -70;
ECa = 50;
EH = -38.8;

%% Somatic currents
INa     = gNa  * m^3 * h * (Vs-ENa);
INaP_s  = gNaP_s * mNaP * hNaP * (Vs-ENa);
IK      = gK   * n^4 * (Vs-EK);
Ileak_s = gleak * (Vs-Eleak);

IT_s = gT_s * mT * hT * (Vs-ECa);
IN_s = gN_s * mN * hN * (Vs-ECa);
IP_s = gP_s * mP * (Vs-ECa);
ICa_s = IT_s + IN_s + IP_s;

IA    = gA   * mA * hA * (Vs-EK);
IH    = gH   * mH * (Vs-EH);
ISK_s = gSK_s * zSK_s^2 * (Vs-EK);

%% Dendritic currents
Ileak_dcur = gleak_d * (Vd-Eleak);

IT_d = gT_d * mT_d * hT_d * (Vd-ECa);
IN_d = gN_d * mN_d * hN_d * (Vd-ECa);
ICa_d = IT_d + IN_d;% + IP_d;

ISK_d = gSK_d * zSK_d^2 * (Vd-EK);

%% Axonal currents
INa_ax   = gNa_ax  * m_ax^3 * h_ax * (Va-ENa);
IK_ax    = gK_ax   * n_ax^4 * (Va-EK);
INaP_ax  = gNaP_ax * mNaP_ax * hNaP_ax * (Va-ENa);
Ileak_ax = gleak_ax * (Va-Eleak);

%% Axial coupling
Iax_sa = g_c_sa * (Vs - Va);   % soma -> axon
Iax_as = g_c_sa * (Va - Vs);   % axon -> soma
Iax_sd = g_c_ds * (Vs - Vd);   % soma -> dendrite
Iax_ds = g_c_ds * (Vd - Vs);   % dendrite -> soma

%% Soma gating
dm = ((1/(1+exp(-(Vs+34)/8.5)))-m)/0.1;

tau_h = 1 + 3.5/(exp((Vs+35)/4) + exp(-(Vs+35)/25));
dh = ((1/(1+exp((Vs+44.1)/7)))-h)/tau_h;

dmNaP = (1/(1+exp(-(Vs+52)/4)) - mNaP)/0.1;
dhNaP = (1/(1+exp((Vs+65)/5)) - hNaP)/tau_hNaP_s;

tau_n = 0.01 + 2.5/(exp((Vs+30)/40) + exp(-(Vs+30)/50));
dn = ((1/(1+exp(-(Vs+30)/25))) - n)/tau_n;

dmT = ((1/(1+exp(-(Vs+38)/5)))-mT)/(2+5/(exp((Vs+28)/25)+exp(-(Vs+28)/70)));
dhT = ((1/(1+exp((Vs+70.1)/7)))-hT)/(1+20/(exp((Vs+70)/65)+exp(-(Vs+70)/65)));

dmN = ((1/(1+exp(-(Vs+30)/6))) - mN)/5;
dhN = ((1/(1+exp((Vs+70)/3))) - hN)/25;

dmP = ((1/(1+exp(-(Vs+17)/3))) - mP)/10;

dmA = ((1/(1+exp(-(Vs+27)/16))) - mA)/(0.37+1/(exp((Vs+40)/5)+exp(-(Vs+74)/7.5)));
dhA = ((1/(1+exp((Vs+80)/11))) - hA)/20;

dmH = ((1/(1+exp(-(Vs+79.8)/5.3))) - mH)/(50+475/(exp((Vs+70)/11)+exp(-(Vs+70)/11)));

dCa_s = K1*ICa_s - K2*Ca_s;
Ca_s_safe = max(Ca_s, 1e-6);
dzSK_s = ((1/(1 + (0.003/Ca_s_safe)^2)) - zSK_s)/1.5;

%% Dendrite gating
dmT_d = ((1/(1+exp(-(Vd+38)/5))) - mT_d) / (2+5/(exp((Vd+28)/25)+exp(-(Vd+28)/70)));
dhT_d = ((1/(1+exp((Vd+70.1)/7))) - hT_d) / (1+20/(exp((Vd+70)/65)+exp(-(Vd+70)/65)));

dmN_d = ((1/(1+exp(-(Vd+30)/6))) - mN_d)/5;
dhN_d = ((1/(1+exp((Vd+70)/3))) - hN_d)/25;

dCa_d = K1*ICa_d - K2*Ca_d;
Ca_d_safe = max(Ca_d, 1e-6);
dzSK_d = ((1/(1 + (0.003/Ca_d_safe)^2)) - zSK_d)/1.5;

%% Axon gating
dm_ax = ((1/(1+exp(-(Va+38)/7)))-m_ax)/0.05;
dh_ax = (1/(1+exp((Va+55)/7)) - h_ax)/1.0;
dn_ax = (1/(1+exp(-(Va+30)/10)) - n_ax)/1.0;

dmNaP_ax = (1/(1+exp(-(Va+52)/4)) - mNaP_ax)/0.1;
dhNaP_ax = (1/(1+exp((Va+56)/2)) - hNaP_ax)/tau_hNaP_ax;

%% Voltage equations
dVs = (-(INa + INaP_s + IK + Ileak_s + ICa_s + ISK_s + IA + IH) ...
    - Iax_sa - Iax_sd + Iext_s ) / Cm_s;

dVa = (-(INa_ax + IK_ax + INaP_ax + Ileak_ax) ...
    - Iax_as ) / Cm_a;

dVd = (-(Ileak_dcur + ICa_d + ISK_d) ...
    - Iax_ds + Iext_d ) / Cm_d;

%% Pack derivatives (29 states)
dydt = [dVd; dVs; dVa; ...
    dm; dh; dmNaP; dhNaP; dn; ...
    dmT; dhT; dmP; dmN; dhN; ...
    dzSK_s; dmA; dhA; dmH; dCa_s; ...
    dmT_d; dhT_d; dmN_d; dhN_d; ...
    dzSK_d; dCa_d; ...
    dmNaP_ax; dhNaP_ax; ...
    dm_ax; dh_ax; dn_ax];

end