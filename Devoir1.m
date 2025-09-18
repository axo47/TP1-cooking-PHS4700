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
    [pcm, partieAvion] = calculCentreMasse(P);
    pcm = pcm - posA;

    % -- direction de poussée (axe x du corps), rotation autour de y --
    ux = [1;0;0];
    if abs(ar) > eps
        R = [ cos(ar), 0,  sin(ar); ...
              0,       1,  0; ...
             -sin(ar), 0,  cos(ar) ];
        for k = 1:numel(partieAvion)
            partieAvion(k).r   = R*(partieAvion(k).r - posA) + posA;
            partieAvion(k).Icm = R * partieAvion(k).Icm * R.';  % I' = R I R^T
        end
        ux = R*ux;
    end

    % -- récup centres de masse utiles --
    names = {partieAvion.name};

    for i = 1:length(names)
        fprintf('  %d: "%s" (classe: %s)\n', i, names{i}, class(names{i}));
    end
    
    % Vérifiez les comparaisons
    comparison = strcmp(names,"moteur_droit");
    fprintf('Résultat strcmp pour moteur_droit: %s\n', mat2str(comparison));
    fprintf('Somme des matches: %d\n', sum(comparison));

    centreMasseMoteurDroit  = partieAvion(strcmp(names,'moteur_droit')).r;
    centreMasseMoteurGauche = partieAvion(strcmp(names,'moteur_gauche')).r;
    centreMasseAileDroite   = partieAvion(strcmp(names,'aile_droite')).r;
    centreMasseAileGauche   = partieAvion(strcmp(names,'aile_gauche')).r;

    % -- points d'application des forces --
    centreMasseMoteurDroitApplique  = centreMasseMoteurDroit  - (P.LONGUEUR_MOTEUR/2)*ux;
    centreMasseMoteurGaucheApplique = centreMasseMoteurGauche - (P.LONGUEUR_MOTEUR/2)*ux;
    centreMasseForcePorteeApplique  = 0.5*(centreMasseAileDroite + centreMasseAileGauche);
    centreMasseForcePorteeApplique  = [centreMasseForcePorteeApplique(1); centreMasseForcePorteeApplique(2); 0];

    pts = struct('pos',{ ...
        centreMasseMoteurDroitApplique, ...
        centreMasseMoteurGaucheApplique, ...
        centreMasseForcePorteeApplique ...
    });

    % -- forces (colonnes) --
    forceMoteurDroit  = Forces(1)*ux;
    forceMoteurGauche = Forces(2)*ux;
    ForcePorte        = [0;0;Forces(3)];
    F = [forceMoteurDroit, forceMoteurGauche, ForcePorte];

    % -- inertie totale au pcm (axes parallèles) --
    MI = calculInertie(partieAvion, pcm);

    % -- dynamique d'Euler : τ = Iα + ω×(Iω) -> α = I^{-1}(τ - ω×(Iω)) --
    momentTotal     = calculMomentForce(pts, F, pcm);
    momentCinetique = MI * va;
    aa = MI \ (momentTotal + cross(momentCinetique, va)); % == τ - ω×(Iω)
end


%% ============================= SOUS-FONCTIONS ============================

function momentTotal = calculMomentForce(vecteurPosForce, F, pcm)
    % Somme des moments τ = Σ (r_i - pcm) × F_i
    momentTotal = zeros(3,1);
    for k = 1:numel(vecteurPosForce)
        momentTotal = momentTotal + cross( (vecteurPosForce(k).pos - pcm), F(:,k) );
    end
end


function [pcm, partieAvion] = calculCentreMasse(P)
    % Centres de masse individuels (dans le repère monde)
    cabineCentreMasse      =  [ P.LONGUEUR_FUSELAGE ; 0;P.EPAISSEUR_AILE+RAYON_CABINE] + (1/4) * [P.HAUTEUR_CABINE;0;0]  ;

    fuselageCentreMasse    =  [P.LONGUEUR_FUSELAGE/2; 0; EPAISSEUR_AILE+RAYON_CABINE];

    ailesCentreMasse = [P.CENTRE_X_AILE; 0; P.EPAISSEUR_AILE/2];

    moteurCentreMasse= [5;0; 2*(P.RAYON_FUSELAGE + P.EPAISSEUR_AILE)];

    AileronCentreMasse     = [P.LARGEUR_AILERON/2; 0; (2*P.RAYON_FUSELAGE)+P.EPAISSEUR_AILE];

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
        struct('name', 'cabine',        'm',P.MASSE_CABINE,   'r',cabineCentreMasse,       'Icm',P.INERTIE_CENTRE_CABINE), ...
        struct('name','fuselage',      'm',P.MASSE_FUSELAGE, 'r',fuselageCentreMasse,     'Icm',P.INTERTIE_CENTRE_FUSELAGE), ...
        struct('name','aile',   'm',P.MASSE_AILE,     'r',ailesCentreMasse,  'Icm',     P.INTERTIE_CENTRE_AILE), ...
        struct('name','moteur',  'm',P.MASSE_MOTEUR,   'r',moteurCentreMasse, 'Icm',P.INTERTIE_CENTRE_MOTEUR), ...
        struct('name','aileron',       'm',P.MASSE_AILERON,  'r',AileronCentreMasse,      'Icm',P.INTERTIE_CENTRE_AILERON) ...
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
