%% ============================ FONCTION PRINCIPALE ========================
% function [pcm, MI, aa] = Dev1(posA, ar, va, Forces)
% posA   : point A (3x1) 
% ar     : angle de rotation autour de y (rad)
% va     : vitesse angulaire ω (3x1)
% Forces : [F_moteurDroit; F_moteurGauche; Portance] (3x1 ou 1x3)
function [pcm, MI, aa] = Devoir1(posA, ar, va, Forces)
    % -- paramètres (constantes + inerties) --
    P = avionParams();  % regroupe toutes tes constantes d'origine

    % -- centre de masse global + sous-parties --
    [pcm_local, partieAvion] = calculCentreMasse(P);

    R = [ cos(ar), 0,  sin(ar); ...
            0,       1,  0; ...
        -sin(ar), 0,  cos(ar) ];

    r_0_global = posA - R * [ P.HAUTEUR_CABINE + P.LONGUEUR_FUSELAGE ; 0 ; P.EPAISSEUR_AILE + P.RAYON_CABINE];
    pcm = r_0_global + R*pcm_local;


    % -- points d'application des forces --
    positionMoteurDroitApplique  =  R * [5 - P.LONGUEUR_MOTEUR/2; -(P.RAYON_FUSELAGE+P.RAYON_MOTEUR); P.RAYON_FUSELAGE+P.EPAISSEUR_AILE];
    positionMoteurGaucheApplique =  R * [5 - P.LONGUEUR_MOTEUR/2; (P.RAYON_FUSELAGE+P.RAYON_MOTEUR); P.RAYON_FUSELAGE+P.EPAISSEUR_AILE];
    positionforcePorteeApplique  =  R * [P.CENTRE_X_AILE; 0; 0]; 


    pcm_moment = R*pcm_local;

    
    pts(1).pos = positionMoteurDroitApplique;
    pts(2).pos = positionMoteurGaucheApplique;
    pts(3).pos = positionforcePorteeApplique;


    % -- forces (colonnes) --
    forceMoteurDroit  = R * Forces(1)*[1;0;0];
    forceMoteurGauche = R * Forces(2)*[1;0;0];
    ForcePorte        = [0;0;Forces(3)];
    F = [forceMoteurDroit, forceMoteurGauche, ForcePorte];

    % -- inertie totale au pcm (axes parallèles) --
    MI = calculInertie(partieAvion, pcm_local, ar);

    % -- dynamique d'Euler : τ = Iα + ω×(Iω) -> α = I^{-1}(τ - ω×(Iω)) --
    momentTotal     = calculMomentForce(pts, F, pcm_moment);
    momentCinetique = MI * va;


    aa = MI \ (momentTotal - cross(va,momentCinetique)); % == τ - ω×(Iω)
end


%% ============================= SOUS-FONCTIONS ============================

function momentTotal = calculMomentForce(vecteurPosForce, F, pcm)
    % Somme des moments τ = Σ (r_i - pcm) × F_i

    momentTotal = zeros(3,1);

    for k = 1:length(vecteurPosForce)

        % Récupération des vecteurs
        r = vecteurPosForce(k).pos;
        force_k = F(:,k);

        % Calcul du bras de levier
        bras_de_levier = r - pcm;


        % Calcul du moment individuel
        moment_k = cross(bras_de_levier, force_k);

        % Ajout au total
        momentTotal = momentTotal + moment_k;
    end

end



function [pcm, partieAvion] = calculCentreMasse(P)
    % Centres de masse individuels (dans le repère monde)
    cabineCentreMasse      =  [ P.LONGUEUR_FUSELAGE ; 0;P.EPAISSEUR_AILE+P.RAYON_CABINE] + (1/4) * [P.HAUTEUR_CABINE;0;0]  ;
    fuselageCentreMasse    =  [P.LONGUEUR_FUSELAGE/2; 0; P.EPAISSEUR_AILE+P.RAYON_CABINE];
    
    ailesCentreMasse = [P.CENTRE_X_AILE; 0; P.EPAISSEUR_AILE/2];
    aileDroiteCentreMasse = [P.CENTRE_X_AILE; -P.CENTRE_Y_AILE ; P.EPAISSEUR_AILE/2];
    aileGaucheCentreMasse = [P.CENTRE_X_AILE; P.CENTRE_Y_AILE ; P.EPAISSEUR_AILE/2];
    
    moteurCentreMasse= [5;0; (P.RAYON_FUSELAGE + P.EPAISSEUR_AILE)];
    moteurGaucheCentreMasse = [5; (P.RAYON_FUSELAGE + P.RAYON_MOTEUR) ; (P.RAYON_FUSELAGE + P.EPAISSEUR_AILE)];
    moteurDroiteCentreMasse = [5; -(P.RAYON_FUSELAGE + P.RAYON_MOTEUR) ; (P.RAYON_FUSELAGE + P.EPAISSEUR_AILE)];


    AileronCentreMasse     = [P.LARGEUR_AILERON/2; 0; (2*P.RAYON_FUSELAGE)+P.EPAISSEUR_AILE + (P.HAUTEUR_AILERON/2)];

    % Barycentre global (ajout du fuselage)
    Sum =  (cabineCentreMasse       * P.MASSE_CABINE) ...
         + (fuselageCentreMasse     * P.MASSE_FUSELAGE) ...
         + (ailesCentreMasse  * 2 * P.MASSE_AILE) ...
         + (moteurCentreMasse * 2 * P.MASSE_MOTEUR) ...
         + (AileronCentreMasse      * P.MASSE_AILERON);

    masseTotale = P.MASSE_CABINE + P.MASSE_FUSELAGE + 2*P.MASSE_AILE + 2*P.MASSE_MOTEUR + P.MASSE_AILERON;
    pcm = Sum / masseTotale;

    % Liste des sous-parties (positions au CM propre + inertie au CM propre)
    partieAvion = [ ...
        struct('name','cabine','m',P.MASSE_CABINE,'r',cabineCentreMasse,'Icm',P.INERTIE_CENTRE_CABINE), ...
        struct('name','fuselage','m',P.MASSE_FUSELAGE,'r',fuselageCentreMasse,'Icm',P.INTERTIE_CENTRE_FUSELAGE), ...        struct('name','ailes','m',P.MASSE_AILE,'r',ailesCentreMasse,'Icm',P.INTERTIE_CENTRE_AILE), ...
        struct('name', 'moteur_droite', 'm', P.MASSE_MOTEUR,'r', moteurDroiteCentreMasse,'Icm', P.INERTIE_CENTRE_MOTEUR), ...
        struct('name', 'moteur_gauche', 'm', P.MASSE_MOTEUR,'r', moteurGaucheCentreMasse,'Icm', P.INERTIE_CENTRE_MOTEUR), ...
        struct('name','aileron','m',P.MASSE_AILERON,'r',AileronCentreMasse,'Icm',P.INERTIE_CENTRE_AILERON), ...
        struct('name', 'aile_droite', 'm', P.MASSE_AILE,'r', aileDroiteCentreMasse, 'Icm', P.INERTIE_CENTRE_AILE), ...
        struct('name', 'aile_gauche', 'm', P.MASSE_AILE,'r', aileGaucheCentreMasse, 'Icm', P.INERTIE_CENTRE_AILE)...
    ];

end


function inertieTotal = calculInertie(partieAvion, pcm, ar)
    % Théorème des axes parallèles sous forme matricielle
    inertieTotal = zeros(3,3);

    for k = 1:numel(partieAvion)
        r = partieAvion(k).r;
        d = pcm - r;           % vecteur de r_k vers pcm
        dx = d(1); dy = d(2); dz = d(3);


        % Matrice T(d) = (||d||^2)I - d d^T
        T = [ dy^2 + dz^2,  -dx*dy,     -dx*dz; ...
              -dx*dy,        dx^2+dz^2, -dy*dz; ...
              -dx*dz,       -dy*dz,      dx^2+dy^2 ];

        inertieTotal = inertieTotal + partieAvion(k).Icm + partieAvion(k).m * T;
    end

    R = [ cos(ar), 0,  sin(ar); ...
            0,       1,  0; ...
        -sin(ar), 0,  cos(ar) ];

    inertieTotal = R * inertieTotal * transpose(R);

end


% ======================= PARAMÈTRES (constantes) ==========================
function P = avionParams()
    % --- Cabine ---
    P.MASSE_CABINE   = 1700;
    P.RAYON_CABINE   = 1.345;
    P.HAUTEUR_CABINE = 3.82;
    P.INERTIE_CENTRE_CABINE = P.MASSE_CABINE * [ ...
         (3*(P.RAYON_CABINE^2))/10, 0, 0; ...
        0, (12*(P.RAYON_CABINE^2) + 3*(P.HAUTEUR_CABINE^2))/80, 0; ...
        0, 0, (12*(P.RAYON_CABINE^2) + 3*(P.HAUTEUR_CABINE^2))/80 ...
    ];

    % --- Fuselage ---
    P.LONGUEUR_FUSELAGE = 22.95;
    P.RAYON_FUSELAGE    = 1.345;
    P.MASSE_FUSELAGE    = 15100;
    P.INTERTIE_CENTRE_FUSELAGE = [ ...
        (P.MASSE_FUSELAGE/2)*(P.RAYON_FUSELAGE^2) ,  0, 0; ...
        0, (P.MASSE_FUSELAGE/4)*(P.RAYON_FUSELAGE^2) + (P.MASSE_FUSELAGE/12)*(P.LONGUEUR_FUSELAGE^2), 0; ...
        0, 0,(P.MASSE_FUSELAGE/4)*(P.RAYON_FUSELAGE^2) + (P.MASSE_FUSELAGE/12)*(P.LONGUEUR_FUSELAGE^2), ...
    ];

        % --- Ailes ---
    P.LONGUEUR_AILE  = 10.6; % x
    P.LARGEUR_AILE   = 1.14; % y
    P.EPAISSEUR_AILE = 0.25; % hauteur z
    P.MASSE_AILE     = 3250;
    P.CENTRE_X_AILE  = 10.54;
    P.CENTRE_Y_AILE  = P.LONGUEUR_AILE/2;

    P.INERTIE_CENTRE_AILE = [ ...
        (P.MASSE_AILE/12)*(P.LONGUEUR_AILE^2 + P.EPAISSEUR_AILE^2), 0, 0; ...
        0, (P.MASSE_AILE/12)*(P.LARGEUR_AILE^2 + P.EPAISSEUR_AILE^2), 0; ...
        0, 0, (P.MASSE_AILE/12)*(P.LARGEUR_AILE^2   + P.LONGUEUR_AILE^2) ...
    ];



    % --- Aileron ---
    P.HAUTEUR_AILERON   = 2.1;
    P.LARGEUR_AILERON   = 1.28;
    P.EPAISSEUR_AILERON = 0.07;
    P.MASSE_AILERON     = 500;
    P.INERTIE_CENTRE_AILERON = [ ...
        (P.MASSE_AILERON/12)*(P.HAUTEUR_AILERON^2 + P.EPAISSEUR_AILERON^2) , 0, 0; ...
        0, (P.MASSE_AILERON/12)*(P.LARGEUR_AILERON^2 + P.HAUTEUR_AILERON^2), 0; ...
        0, 0, (P.MASSE_AILERON/12)*(P.LARGEUR_AILERON^2 + P.EPAISSEUR_AILERON^2) ...
    ];

    % --- Moteur ---
    P.LONGUEUR_MOTEUR = 3.68;
    P.RAYON_MOTEUR    = 0.724;
    P.MASSE_MOTEUR    = 1700;
    P.INERTIE_CENTRE_MOTEUR = [ ...
        (P.MASSE_MOTEUR/2)*(P.RAYON_MOTEUR^2), 0, 0; ...
        0, (P.MASSE_MOTEUR/4)*(P.RAYON_MOTEUR^2) + (P.MASSE_MOTEUR/12)*(P.LONGUEUR_MOTEUR^2), 0; ...
        0, 0,  (P.MASSE_MOTEUR/4)*(P.RAYON_MOTEUR^2) + (P.MASSE_MOTEUR/12)*(P.LONGUEUR_MOTEUR^2) ...
    ];

    % positions-type si tu veux les réutiliser ailleurs
    P.POSITION_FORCE_MOTEUR_DROIT = [ 5 - (P.LONGUEUR_MOTEUR/2); +(P.RAYON_FUSELAGE+P.RAYON_MOTEUR); P.RAYON_FUSELAGE + P.EPAISSEUR_AILE ];
    P.POSITION_FORCE_MOTEUR_GAUCHE= [ 5 - (P.LONGUEUR_MOTEUR/2); -(P.RAYON_FUSELAGE+P.RAYON_MOTEUR); P.RAYON_FUSELAGE + P.EPAISSEUR_AILE ];
    P.POSITION_FORCE_PORTE        = [ P.CENTRE_X_AILE; P.CENTRE_Y_AILE; 0 ];
end