function FEM = Cal_FGM_HSDT(FEM)
global h 

% ===== Material matrix =====
Ab=zeros(3,3); Bb=Ab; Db=Ab; Eb=Ab; Fb=Ab; Hb=Ab;
As=zeros(2,2); Bs=zeros(2,2); Ds=zeros(2,2); I=zeros(1,6);

% ===== Material properties ======
if FEM.Material == 0 % Isotropic material
    Em = 30e6; Num = 0.3; Rhom = 2300;  %
    Ec = 30e6; Nuc = 0.3; Rhoc = 2300; %
    FEM.Em = Em; FEM.Num = Num;  FEM.Rho = Rhom;
elseif FEM.Material == 1
    Em = 70e9; Num = 0.3; Rhom = 2707;  % Al
    Ec = 380e9; Nuc = 0.3; Rhoc = 3800; % Al2O3
    FEM.Em = Em; FEM.Num = Num; FEM.Rhom = Rhom;
    FEM.Ec = Ec; FEM.Nuc = Nuc; FEM.Rhoc = Rhoc;
elseif FEM.Material == 2 % Al/ZrO2
    Em = 70e9; Num = 0.3; Rhom = 2707;  %
    Ec = 151e9; Nuc = 0.3; Rhoc = 5700; %   ZrO2-2
    FEM.Em = Em; FEM.Num = Num; FEM.Rhom = Rhom;
    FEM.Ec = Ec; FEM.Nuc = Nuc; FEM.Rhoc = Rhoc;
elseif FEM.Material == 3 % SUS304/Si3N4
    Em = 201.04e9; Num = 0.3;    Rhom = 8166;  % SUS304
    Ec = 348.43e9; Nuc = 0.3;    Rhoc = 2370;  % Si3N4
    FEM.Em = Em; FEM.Num = Num; FEM.Rhom = Rhom;
    FEM.Ec = Ec; FEM.Nuc = Nuc; FEM.Rhoc = Rhoc;
end
% ==== Get integration points ====
x1 = 1; x2 = -1 ;
[Qg,Wg] = gauleg(x1,x2,30) ;
% ==== Core layer ====
zk = [-h/2  h/2];
n = FEM.IndexMat;

for igp = 1 : size(Wg,1) % z-direction
    pt = Qg(igp,:) ;
    gpt = (zk(2) + zk(1))/2 + (zk(2) - zk(1))/2*pt;  % value of z-direction
    wt = (zk(2) - zk(1))/2*Wg(igp) ;
    
    Vc =  (0.5 + gpt/h)^n;
    
    if FEM.Law == 0 % Pover law
            E_z = Em + Vc*(Ec - Em);
            Nu_z = Num + Vc*(Nuc - Num);
            Rho_z = Rhom + Vc*(Rhoc - Rhom);
    elseif FEM.Law == 1 % Mori_Tanaka
        Vm = 1 - Vc;
        Gc = Ec/(2*(1 + Nuc)) ;      Gm = Em/(2*(1 + Num)) ;           % shear modulus
        shearc = Ec/3/(1 - 2*Nuc);   shearm = Em/3/(1 - 2*Num);        % bulk modulus
        f1 =  Gm*(9*shearm + 8*Gm)/(6*(shearm + 2*Gm)) ;
        Geff = Gm + (Gc - Gm)*Vc/(1 + Vm*(Gc - Gm)/(Gm + f1)) ;
        shearcm = shearc - shearm ;
        Sheareff = shearm + shearcm*Vc/(1 + Vm*shearcm/(shearm + (4/3)*Gm)) ;
        % Effective Young moduls and Poisson ratio and density
        E_z = 9*Sheareff*Geff/(3*Sheareff + Geff) ;
        Nu_z = (3*Sheareff - 2*Geff)/(2*(3*Sheareff + Geff)) ;
        Rhocm = Rhoc - Rhom ;
        Rho_z = Rhom + Rhocm*Vc ;          % effective density
    end
    
        % ==== Bending ====
        Cb = E_z/(1-Nu_z^2)*[1   Nu_z   0;
                            Nu_z  1   0;
                            0   0  (1- Nu_z)/2] ;
        % ==== Shear ====
        Cs = E_z/2/(1+ Nu_z)*[1 0;
                             0 1];

        [ff,dff] = Dist_Func(h,gpt,FEM.DisFunc);

        Ab = Ab + Cb*wt ;       Bb = Bb + Cb*gpt*wt ;    Db = Db + Cb*gpt^2*wt ;
        Eb = Eb + Cb*ff*wt ;    Fb = Fb + Cb*ff*gpt*wt ; Hb = Hb + Cb*ff^2*wt ;

        As = As + Cs*wt;        Bs = Bs + Cs*dff*wt;
        Ds = Ds + Cs*dff^2*wt;
    % ================= inertia terms matrix =======================
    I = I + (Rho_z*wt).* [1 gpt gpt^2 ff ff*gpt ff^2];
end

FEM.DDb = [Ab Bb Eb;Bb Db Fb;Eb Fb Hb];
FEM.DDs = [As Bs; Bs Ds];
I1 = I(1)*eye(3); I2 = I(2)*eye(3); I3 = I(3)*eye(3);
I4 = I(4)*eye(3); I5 = I(5)*eye(3); I6 = I(6)*eye(3);

FEM.Im = [I1 I2 I4; I2 I3 I5;I4 I5 I6];
