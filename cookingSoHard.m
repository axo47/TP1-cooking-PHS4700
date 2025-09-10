%% ============================ FONCTION PRINCIPALE ========================
% function [pcm, MI, aa] = Dev1(posA, ar, va, Forces)
% posA   : point A (3x1) à partir duquel on référence le fuselage/cabine
% ar     : angle de rotation autour de y (rad)
% va     : vitesse angulaire ω (3x1)
% Forces : [F_moteurDroit; F_moteurGauche; Portance] (3x1 ou 1x3)
function [pcm, MI, aa] = Dev1(posA, ar, va, Forces)
    % Rendre les constantes visibles ici (fidèle à ta structure)
    global LONGUEUR_MOTEUR RAYON_FUSELAGE RAYON_MOTEUR EPAISSEUR_AILE

    % Petites sécurités sur les entrées (sans changer ton API)
    assert(isvector(posA) && numel(posA)==3, 'posA doit être 3x1.');
    assert(isvector(va)   && numel(va)==3,   'va (ω) doit être 3x1.');
    Forces = Forces(:);
    assert(numel(Forces)==3, 'Forces doit contenir [F_droit, F_gauche, Portance].');

    % Centre de masse global + liste des sous-parties
    [pcm, partieAvion] = calculCentreMasse(posA);

    % Axe x du repère corps (pour pousser les moteurs)
    ux = [1;0;0];

    % Rotation autour de y si ar ~= 0
    if abs(ar) > eps
        R = [ cos(ar), 0,  sin(ar); ...
              0,       1,  0; ...
             -sin(ar), 0,  cos(ar) ];
        % Rotation positions et inerties autour de pcm
        for k = 1:numel(partieAvion)
            partieAvion(k).r   = pcm + R*(partieAvion(k).r - pcm);
            partieAvion(k).Icm = R * partieAvion(k).Icm * R.';  % I' = R I R^T
        end
        ux = R*ux;  % direction de poussée après rotation
    end

    % Récup des centres de masse utiles
    names = {partieAvion.name};
    centreMasseMoteurDroit = partieAvion(strcmp(names,"moteur_droit")).r;
    centreMasseMoteurGauche= partieAvion(strcmp(names,"moteur_gauche")).r;
    centreMasseAileDroite  = partieAvion(strcmp(names,"aile_droite")).r;
    centreMasseAileGauche  = partieAvion(strcmp(names,"aile_gauche")).r;

    % Points d’application des forces
    centreMasseMoteurDroitApplique  = centreMasseMoteurDroit  - (LONGUEUR_MOTEUR/2)*ux;
    centreMasseMoteurGaucheApplique = centreMasseMoteurGauche - (LONGUEUR_MOTEUR/2)*ux;
    centreMasseForcePorteeApplique  = 0.5*(centreMasseAileDroite + centreMasseAileGauche);
    centreMasseForcePorteeApplique  = [centreMasseForcePorteeApplique(1); centreMasseForcePorteeApplique(2); 0];

    pts = struct('pos',{ ...
        centreMasseMoteurDroitApplique, ...
        centreMasseMoteurGaucheApplique, ...
        centreMasseForcePorteeApplique ...
    });

    % Forces (colonnes)
    forceMoteurDroit  = Forces(1)*ux;
    forceMoteurGauche = Forces(2)*ux;
    ForcePorte        = [0;0;Forces(3)];
    F = [forceMoteurDroit, forceMoteurGauche, ForcePorte];

    % Inertie totale au pcm (théorème des axes parallèles)
    MI = calculInertie(partieAvion, pcm);

    % Moments et accélération angulaire (Euler : τ = Iα + ω×(Iω))
    momentTotal     = calculMomentForce(pts, F, pcm);
    momentCinetique = MI * va;                     % Iω
    aa = MI \ (momentTotal + cross(momentCinetique, va)); 
end


%% ============================= SOUS-FONCTIONS ============================

function momentTotal = calculMomentForce(vecteurPosForce, F, pcm)
    % Somme des moments τ = Σ (r_i - pcm) × F_i
    momentTotal = zeros(3,1);
    for k = 1:numel(vecteurPosForce)
        momentTotal = momentTotal + cross( (vecteurPosForce(k).pos - pcm), F(:,k) );
    end
end


function [pcm, partieAvion] = calculCentreMasse(posA)
    % Rendre visibles toutes les constantes et inerties
    global MASSE_CABINE RAYON_FUSELAGE HAUTEUR_CABINE ...
           LONGUEUR_FUSELAGE MASSE_FUSELAGE ...
           LONGUEUR_AILE LARGEUR_AILE EPAISSEUR_AILE MASSE_AILE CENTRE_X_AILE CENTRE_Y_AILE ...
           LONGUEUR_MOTEUR RAYON_MOTEUR MASSE_MOTEUR ...
           HAUTEUR_AILERON LARGEUR_AILERON EPAISSEUR_AILERON MASSE_AILERON ...
           INERTIE_CENTRE_CABINE INTERTIE_CENTRE_FUSELAGE INTERTIE_CENTRE_AILE INTERTIE_CENTRE_MOTEUR INTERTIE_CENTRE_AILERON

    % Centres de masse individuels (dans le repère monde)
    baseCabine            = [-(3/4)*HAUTEUR_CABINE; 0; 0];
    cabineCentreMasse     = posA + baseCabine;

    fuselageCentreMasse   = posA + [-(HAUTEUR_CABINE + LONGUEUR_FUSELAGE/2); 0; 0];

    ailesCentreMasseDroite= [CENTRE_X_AILE; +CENTRE_Y_AILE; EPAISSEUR_AILE];
    ailesCentreMasseGauche= [CENTRE_X_AILE; -CENTRE_Y_AILE; EPAISSEUR_AILE];

    moteurCentreMasseDroite = [5; +(RAYON_FUSELAGE+RAYON_MOTEUR); RAYON_FUSELAGE + EPAISSEUR_AILE];
    moteurCentreMasseGauche = [5; -(RAYON_FUSELAGE+RAYON_MOTEUR); RAYON_FUSELAGE + EPAISSEUR_AILE];

    AileronCentreMasse    = [LARGEUR_AILERON/2; 0; (2*RAYON_FUSELAGE)+EPAISSEUR_AILERON];

    % Barycentre global
    Sum = (cabineCentreMasse      * MASSE_CABINE) ...
        + (ailesCentreMasseDroite * MASSE_AILE) ...
        + (ailesCentreMasseGauche * MASSE_AILE) ...
        + (moteurCentreMasseDroite* MASSE_MOTEUR) ...
        + (moteurCentreMasseGauche* MASSE_MOTEUR) ...
        + (AileronCentreMasse     * MASSE_AILERON);

    masseTotale = MASSE_CABINE + 2*MASSE_AILE + 2*MASSE_MOTEUR + MASSE_AILERON;
    pcm = Sum / masseTotale;

    % Liste des sous-parties (positions au CM propre + inertie au CM propre)
    partieAvion = [ ...
        struct('name',"cabine",        'm',MASSE_CABINE,   'r',cabineCentreMasse,       'Icm',INERTIE_CENTRE_CABINE), ...
        struct('name',"fuselage",      'm',MASSE_FUSELAGE, 'r',fuselageCentreMasse,     'Icm',INTERTIE_CENTRE_FUSELAGE), ...
        struct('name',"aile_droite",   'm',MASSE_AILE,     'r',ailesCentreMasseDroite,  'Icm',INTERTIE_CENTRE_AILE), ...
        struct('name',"aile_gauche",   'm',MASSE_AILE,     'r',ailesCentreMasseGauche,  'Icm',INTERTIE_CENTRE_AILE), ...
        struct('name',"moteur_droit",  'm',MASSE_MOTEUR,   'r',moteurCentreMasseDroite, 'Icm',INTERTIE_CENTRE_MOTEUR), ...
        struct('name',"moteur_gauche", 'm',MASSE_MOTEUR,   'r',moteurCentreMasseGauche, 'Icm',INTERTIE_CENTRE_MOTEUR), ...
        struct('name',"aileron",       'm',MASSE_AILERON,  'r',AileronCentreMasse,      'Icm',INTERTIE_CENTRE_AILERON) ...
    ];
end


function inertieTotal = calculInertie(partieAvion, pcm)
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
end


% ======================= PARAMÈTRES (constantes) ==========================
function P = avionParams()
    % --- Cabine ---
    P.MASSE_CABINE   = 1700;
    P.RAYON_CABINE   = 1.345;
    P.HAUTEUR_CABINE = 3.82;
    P.INERTIE_CENTRE_CABINE = P.MASSE_CABINE * [ ...
        (12*(P.RAYON_CABINE^2) + 3*(P.HAUTEUR_CABINE^2))/80, 0, 0; ...
        0, (12*(P.RAYON_CABINE^2) + 3*(P.HAUTEUR_CABINE^2))/80, 0; ...
        0, 0, (3*(P.RAYON_CABINE^2))/10 ...
    ];

    % --- Fuselage ---
    P.LONGUEUR_FUSELAGE = 22.95;
    P.RAYON_FUSELAGE    = 1.345;
    P.MASSE_FUSELAGE    = 15100;
    P.INTERTIE_CENTRE_FUSELAGE = [ ...
        (P.MASSE_FUSELAGE/4)*(P.RAYON_FUSELAGE^2) + (P.MASSE_FUSELAGE/12)*(P.LONGUEUR_FUSELAGE^2), 0, 0; ...
        0, (P.MASSE_FUSELAGE/4)*(P.RAYON_FUSELAGE^2) + (P.MASSE_FUSELAGE/12)*(P.LONGUEUR_FUSELAGE^2), 0; ...
        0, 0, (P.MASSE_FUSELAGE/2)*(P.RAYON_FUSELAGE^2) ...
    ];

    % --- Ailes ---
    P.LONGUEUR_AILE  = 10.6;
    P.LARGEUR_AILE   = 1.14;
    P.EPAISSEUR_AILE = 0.25;
    P.MASSE_AILE     = 3250;
    P.CENTRE_X_AILE  = 10.54;
    P.CENTRE_Y_AILE  = P.RAYON_FUSELAGE + (P.LONGUEUR_AILE/2);
    P.INTERTIE_CENTRE_AILE = [ ...
        (P.MASSE_AILE/12)*(P.LARGEUR_AILE^2   + P.EPAISSEUR_AILE^2), 0, 0; ...
        0, (P.MASSE_AILE/12)*(P.LONGUEUR_AILE^2 + P.EPAISSEUR_AILE^2), 0; ...
        0, 0, (P.MASSE_AILE/12)*(P.LONGUEUR_AILE^2 + P.LARGEUR_AILE^2) ...
    ];

    % --- Aileron ---
    P.HAUTEUR_AILERON   = 2.1;
    P.LARGEUR_AILERON   = 1.28;
    P.EPAISSEUR_AILERON = 0.07;
    P.MASSE_AILERON     = 500;
    P.INTERTIE_CENTRE_AILERON = [ ...
        (P.MASSE_AILERON/12)*(P.LARGEUR_AILERON^2 + P.EPAISSEUR_AILERON^2), 0, 0; ...
        0, (P.MASSE_AILERON/12)*(P.HAUTEUR_AILERON^2 + P.EPAISSEUR_AILERON^2), 0; ...
        0, 0, (P.MASSE_AILERON/12)*(P.HAUTEUR_AILERON^2 + P.LARGEUR_AILERON^2) ...
    ];

    % --- Moteur ---
    P.LONGUEUR_MOTEUR = 3.68;
    P.RAYON_MOTEUR    = 0.724;
    P.MASSE_MOTEUR    = 1700;
    P.INTERTIE_CENTRE_MOTEUR = [ ...
        (P.MASSE_MOTEUR/4)*(P.RAYON_MOTEUR^2) + (P.MASSE_MOTEUR/12)*(P.LONGUEUR_MOTEUR^2), 0, 0; ...
        0, (P.MASSE_MOTEUR/4)*(P.RAYON_MOTEUR^2) + (P.MASSE_MOTEUR/12)*(P.LONGUEUR_MOTEUR^2), 0; ...
        0, 0, (P.MASSE_MOTEUR/2)*(P.RAYON_MOTEUR^2) ...
    ];

    % positions-type si tu veux les réutiliser ailleurs
    P.POSITION_FORCE_MOTEUR_DROIT = [ 5 - (P.LONGUEUR_MOTEUR/2); +(P.RAYON_FUSELAGE+P.RAYON_MOTEUR); P.RAYON_FUSELAGE + P.EPAISSEUR_AILE ];
    P.POSITION_FORCE_MOTEUR_GAUCHE= [ 5 - (P.LONGUEUR_MOTEUR/2); -(P.RAYON_FUSELAGE+P.RAYON_MOTEUR); P.RAYON_FUSELAGE + P.EPAISSEUR_AILE ];
    P.POSITION_FORCE_PORTE        = [ P.CENTRE_X_AILE; P.CENTRE_Y_AILE; 0 ];
end
