function FEM = Cal_TPMS_GPL_E_control(FEM)
%%
global E_Dis e_0 h
global GPL_Dis l_gpl w_gpl t_gpl E_gpl nu_gpl rho_gpl lambda_gpl
global Porous_type Mat_Model

% ===== Material matrix =====
Ab=zeros(3,3); Bb=Ab; Db=Ab; Eb=Ab; Fb=Ab; Hb=Ab;
As=zeros(2,2); Bs=zeros(2,2); Ds=zeros(2,2); I=zeros(1,6);

if Mat_Model == 1
    Material = [200e9 0.3 8000]; % Steel
elseif Mat_Model == 2
    Material = [130e9 0.34 8960]; % Copper
end

% ===== material property======
x1 = 1; x2 = -1 ;
[Qg,Wg] = gauleg(x1,x2,20) ;
zu = h/2;              % upper
zl = -h/2;             % lower
% ===== Young's modulus, poisson's ratio
E_m = Material(1); FEM.Em = E_m;
nu_m = Material(2); FEM.Num = nu_m;
G_m = E_m/(2*(1+nu_m));
rho_m = Material(3);FEM.Rhom = rho_m;

%% Define parameter of TPMS
switch Porous_type
    case {1, 2, 3}
        if Porous_type == 1 % Primitive
            k_e = 0.25; C_1e = 0.317; n_1e = 1.264; n_2e = 2.006;
            k_g = 0.25; C_1g = 0.705; n_1g = 1.189; n_2g = 1.715;
            k_nu = 0.55; a_1 = 0.314; b_1 = -1.004; a_2 = 0.152;
        elseif Porous_type == 2 % Gyroid
            k_e = 0.45; C_1e = 0.596; n_1e = 1.467; n_2e = 2.351;
            k_g = 0.45; C_1g = 0.777; n_1g = 1.544; n_2g = 1.982;
            k_nu = 0.50; a_1 = 0.192; b_1 = -1.349; a_2 = 0.402;
        elseif Porous_type == 3 % IWP
            k_e = 0.35; C_1e = 0.597; n_1e = 1.225; n_2e = 1.782;
            k_g = 0.35; C_1g = 0.529; n_1g = 1.287; n_2g = 2.188;
            k_nu = 0.13; a_1 = 2.597; b_1 = -0.157; a_2 = 0.201;
        end
        C_2e = (C_1e*k_e^(n_1e) - 1)/(k_e^(n_2e) - 1); C_3e = 1 - C_2e;
        C_2g = (C_1g*k_g^(n_1g) - 1)/(k_g^(n_2g) - 1); C_3g = 1 - C_2g;
        d_1 = 0.3 - a_1*exp(b_1*k_nu); b_2 = - a_2*(k_nu + 1); d_2 = 0.3 - a_2*(1)^2 - b_2(1);
end

%% Calculate psi_{e} of porous ditribution 3
M = 0;
if E_Dis == 3
    for igp = 1:size(Wg,1)
        pt = Qg(igp,:) ;

        % map the point to global
        gpt = (zu+zl)/2 + (zu-zl)/2*pt ;                 % value of z-direction
        wt = (zu-zl)/2*Wg(igp) ;

        psi_z = cos(pi*gpt/h);

        % 1: Primitive, 2: Gyroid, 3: IWP, 4: Closed-cell, 5: Open-cell
        switch Porous_type
            case {1, 2, 3}
                if psi_z >= (1 - C_1e*k_e^(n_1e))/e_0
                    e_m = 1/psi_z - 1/psi_z * ((1 - e_0*psi_z)/(C_1e))^(1/n_1e);
                else
                    e_m = 1/psi_z - 1/psi_z * ((1 - e_0*psi_z - C_3e)/(C_2e))^(1/n_2e);
                end
            case 4
                e_m = 1.121/psi_z - 1.121/psi_z * (1 -e_0*psi_z)^(1/2.3);
            case 5
                e_m = 1/psi_z - 1/psi_z * (1 -e_0*psi_z)^(1/2);
        end

        M = M + (1 - e_m*psi_z)*wt;
    end

    switch Porous_type
        case {1, 2, 3}
            psi_e = 1/e_0 - C_1e/e_0 * (M/h)^(n_1e);
            if psi_e < (1 - C_1e*k_e^(n_1e))/e_0
                psi_e = (1-C_3e)/e_0 - C_2e/e_0 * (M/h)^(n_2e);
            end
            if psi_e >= (1 - C_1e*k_e^(n_1e))/e_0
                disp('Error in \psi_e');
                pause
            end
        case 4
            psi_e = 1/e_0 - 1/e_0 * ((M/h + 0.121)/1.121)^(2.3);
        case 5
            psi_e = 1/e_0 - 1/e_0 * (M/h)^(2);
    end
end

%% Calculate Si of GPL distributions
Vtotal_GPL = lambda_gpl*rho_m / (lambda_gpl*rho_m + rho_gpl - lambda_gpl*rho_gpl);
A1 = 0; A2 = 0;

for igp = 1:size(Wg,1)
    pt = Qg(igp,:) ;

    % map the point to global
    gpt = (zu+zl)/2 + (zu-zl)/2*pt ;                 % value of z-direction
    wt = (zu-zl)/2*Wg(igp) ;

    switch E_Dis
        case 1
            psi_z = cos(pi*gpt/h);
        case 2
            psi_z = cos(pi*gpt/2/h + pi/4);
        case 3
            psi_z = psi_e;
    end

    % 1: Primitive, 2: Gyroid, 3: IWP, 4: Closed-cell, 5: Open-cell
    switch Porous_type
        case {1, 2, 3}
            if psi_z >= (1 - C_1e*k_e^(n_1e))/e_0
                e_m = 1/psi_z - 1/psi_z * ((1 - e_0*psi_z)/(C_1e))^(1/n_1e);
            else
                e_m = 1/psi_z - 1/psi_z * ((1 - e_0*psi_z - C_3e)/(C_2e))^(1/n_2e);
            end
        case 4
            e_m = 1.121/psi_z - 1.121/psi_z * (1 -e_0*psi_z)^(1/2.3);
        case 5
            e_m = 1/psi_z - 1/psi_z * (1 -e_0*psi_z)^(1/2);
    end

    switch GPL_Dis
        case 1
            Vgpl_z = (1 - cos(pi*gpt/h));
        case 2
            Vgpl_z = (1 - cos(pi*gpt/2/h + pi/4));
        case 3
            Vgpl_z = 1;
    end

    A1 = A1 + (1 - e_m*psi_z)*wt;
    A2 = A2 + (Vgpl_z)*(1 - e_m*psi_z)*wt;
end
S = Vtotal_GPL*A1/A2;

%% Calculate material matrices
xi_l = 2*l_gpl/t_gpl; xi_w = 2*w_gpl/t_gpl;
eta_l = ((E_gpl/E_m) - 1)/((E_gpl/E_m) + xi_l); eta_w = ((E_gpl/E_m) - 1)/((E_gpl/E_m) + xi_w);

for igp = 1:size(Wg,1)
    pt = Qg(igp,:) ;

    % map the point to global
    gpt = (zu+zl)/2 + (zu-zl)/2*pt ;                 % value of z-direction
    wt = (zu-zl)/2*Wg(igp) ;

    % GPLR base material
    switch GPL_Dis
        case 1
            Vgpl_z = S * (1 - cos(pi*gpt/h));
        case 2
            Vgpl_z = S * (1 - cos(pi*gpt/2/h + pi/4));
        case 3
            Vgpl_z = S;
    end

    V_m = 1 - Vgpl_z;
    E_1 = 3/8*(1 + xi_l*eta_l*Vgpl_z)/(1-eta_l*Vgpl_z)*E_m + 5/8*(1 + xi_w*eta_w*Vgpl_z)/(1-eta_w*Vgpl_z)*E_m;
    nu_1 = nu_gpl*Vgpl_z + nu_m*V_m;
    G_1 = E_1 / (2*(1+nu_1));
    rho_1 = rho_gpl*Vgpl_z + rho_m*V_m;

    % Including the porosity
    switch E_Dis
        case 1
            psi_z = cos(pi*gpt/h);
        case 2
            psi_z = cos(pi*gpt/2/h + pi/4);
        case 3
            psi_z = psi_e;
    end

    % 1: Primitive, 2: Gyroid, 3: IWP, 4: Closed-cell, 5: Open-cell
    E_z = E_1 * (1 - e_0*psi_z);
    switch Porous_type
        case {1, 2, 3}
            if psi_z >= (1 - C_1e*k_e^(n_1e))/e_0
                e_m = 1/psi_z - 1/psi_z * ((1 - e_0*psi_z)/(C_1e))^(1/n_1e);
            else
                e_m = 1/psi_z - 1/psi_z * ((1 - e_0*psi_z - C_3e)/(C_2e))^(1/n_2e);
            end

            RD_z = 1 - e_m*psi_z;

            e_z = (RD_z <= k_e) * (C_1e*RD_z^(n_1e)) + ...
                (RD_z >  k_e) * (C_2e*RD_z.^(n_2e) + C_3e);
            g_z = (RD_z <= k_g) * (C_1g*RD_z^(n_1g)) + ...
                (RD_z >  k_g) * (C_2g*RD_z.^(n_2g) + C_3g);
            nu_z =(RD_z <= k_nu) * (a_1.*exp(b_1*RD_z) + d_1) + ...
                (RD_z >  k_nu) * (a_2*RD_z^(2) + b_2*RD_z + d_2);

            %             E_z = E_1*e_z;
            G_z = G_1*g_z;
            %             rho_z = rho*RD;
            if abs(E_z - E_1*e_z) > 1e-3 % Check code
                disp(RD_z)
            end
        case 4
            e_m = 1.121/psi_z - 1.121/psi_z * (1 -e_0*psi_z)^(1/2.3);
            p_star = 1.121*(1 - (1 -e_0*psi_z)^(1/2.3));
            nu_z = 0.221*p_star + nu_1*(0.342*p_star^(2) - 1.21*p_star + 1);
            G_z = E_z / (2*(1+nu_z));
        case 5
            e_m = 1/psi_z - 1/psi_z * (1 -e_0*psi_z)^(1/2);
            nu_z = 1/3;
            G_z = E_z /(2*(1+nu_z));
            %             G_z = G_1 * (1-e_0*psi_z);
            %             if abs(G_z - G_1 * (1-e_0*psi_z)) > 1e-3 % Check code
            %                 disp(G_z - G_1 * (1-e_0*psi_z))
            %             end
    end
    rho_z = rho_1 * (1 - e_m*psi_z);

    % % ============= Bending-extension stiffness matrix ===================
    Qb = [E_z/(1-nu_z^2) E_z*nu_z/(1-nu_z^2) 0;
        E_z*nu_z/(1-nu_z^2) E_z/(1-nu_z^2) 0;
        0 0 G_z];
    % shear stiffness
    Qs =G_z*[1 0;
        0 1];
    
    % choos model for displacement field
     [gg,dgg] = Dist_Func(h,gpt,FEM.DisFunc);

    Ab = Ab + Qb*wt ;        Bb = Bb + Qb*gpt*wt ;     Db = Db + Qb*gpt^2*wt ;
    Eb = Eb + Qb*gg*wt ;     Fb = Fb + Qb*gg*gpt*wt ;  Hb = Hb + Qb*gg^2*wt ;

    As = As + Qs*wt;        Bs = Bs + Qs*dgg*wt;
    Ds = Ds + Qs*dgg^2*wt;
    % ================= Inertia terms matrix =======================
    I = I + (rho_z*wt).* [1 gpt gpt^2 gg gg*gpt gg^2];
end

FEM.DDb = [Ab Bb Eb;Bb Db Fb;Eb Fb Hb];
FEM.DDs = [As Bs; Bs Ds];
I1 = I(1)*eye(3); I2 = I(2)*eye(3); I3 = I(3)*eye(3);
I4 = I(4)*eye(3); I5 = I(5)*eye(3); I6 = I(6)*eye(3);

FEM.Im = [I1 I2 I4; I2 I3 I5;I4 I5 I6];


