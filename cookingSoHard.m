MASSE_CABINE = 1700
RAYON_CABINE = 1.345
HAUTEUR_CABINE = 3.82
INERTIE_CENTRE_CABINE = MASSE_CABINE * [ (12*(RAYON_CABINE^2) + (3*(HAUTEUR_CABINE^2)))/80 , 0 , 0;
                                0 , 12*(RAYON_CABINE^2) + (3*(HAUTEUR_CABINE^2))/80 , 0 ;
                                0 , 0 , (3*(RAYON_CABINE^2))/10 ]


LONGUEUR_FUSELAGE = 22.95
RAYON_FUSELAGE = 1.345
MASSE_FUSELAGE = 15100
INTERTIE_CENTRE_FUSELAGE = [ ((MASSE_FUSELAGE)/4) * (RAYON_FUSELAGE^2) + ((MASSE_FUSELAGE)/12 * (LONGUEUR_FUSELAGE^2)) , 0 , 0;
                                0 , (((MASSE_FUSELAGE)/4) * RAYON_FUSELAGE^2) + ((MASSE_FUSELAGE)/12 * (LONGUEUR_FUSELAGE^2))  , 0 ;
                                0 , 0 , (MASSE_FUSELAGE/2)*(RAYON_FUSELAGE^2) ] 

LONGUEUR_AILE = 10.6
LARGEUR_AILE = 1.14
EPAISSEUR_AILE = 0.25
MASSE_AILE = 3250
CENTRE_X_AILE = 10.54
CENTRE_Y_AILE = RAYON_FUSELAGE+(LONGUEUR_AILE/2)
INTERTIE_CENTRE_AILE = [ (((MASSE_AILE)/12) *  ((LARGEUR_AILE^2) + (EPAISSEUR_AILE^2)))  , 0 , 0;
                                0 , (((MASSE_AILE)/12) *  ((LONGUEUR_AILE^2) + (EPAISSEUR_AILE^2)))  , 0 ;
                                0 , 0 , (((MASSE_AILE)/12) *  ((LONGUEUR_AILE^2) + (LARGEUR_AILE^2))) 
                ]


HAUTEUR_AILERON = 2.1
LARGEUR_AILERON = 1.28
EPAISSEUR_AILERON = 0.07
MASSE_AILERON = 500
INTERTIE_CENTRE_AILERON = [ (((MASSE_AILERON)/12) *  ((LARGEUR_AILERON^2) + (EPAISSEUR_AILERON^2)))  , 0 , 0;
                                0 , (((MASSE_AILERON)/12) *  ((HAUTEUR_AILERON^2) + (EPAISSEUR_AILERON^2)))  , 0 ;
                                0 , 0 , (((MASSE_AILERON)/12) *  ((HAUTEUR_AILERON^2) + (LARGEUR_AILERON^2))) 
                          ]



LONGUEUR_MOTEUR = 3.68
RAYON_MOTEUR = 0.724
MASSE_MOTEUR = 1700
INTERTIE_CENTRE_MOTEUR =  [ (((MASSE_MOTEUR)/4) * RAYON_MOTEUR^2) + ((MASSE_MOTEUR)/12 * (LONGUEUR_MOTEUR^2)) , 0 , 0;
                                0 , (((MASSE_MOTEUR)/4) * RAYON_MOTEUR^2) + ((MASSE_MOTEUR)/12 * (LONGUEUR_MOTEUR^2))  , 0 ;
                                0 , 0 , (MASSE_MOTEUR/2)*(RAYON_MOTEUR^2) ] 



POSITION_FORCE_MOTEUR_DROIT = [5 - (LONGUEUR_MOTEUR/2) ; (RAYON_FUSELAGE+RAYON_MOTEUR); RAYON_FUSELAGE + EPAISSEUR_AILE]
POSITION_FORCE_MOTEUR_GAUCHE = [5 - (LONGUEUR_MOTEUR/2) ; -(RAYON_FUSELAGE+RAYON_MOTEUR); RAYON_FUSELAGE + EPAISSEUR_AILE]
POSITION_FORCE_PORTE = [CENTRE_X_AILE ;  CENTRE_Y_AILE;  0]





function [pcm MI aa] =  Dev1(posA,ar,va,Forces) 
[pcm,partieAvion] = calculCentreMasse(posA);


ux  = [1;0;0];
if abs(ar) ~= 0
    R = [    cos(ar), 0, sin(ar);
             0, 1,0;
             -sin(ar), 0 , cos(ar)   
    ]

    for k = 1:numel(partieAvion)
        partieAvion(k).r = pcm + R * (partieAvion(k).r-pcm);
        partieAvion(k).Icm = R * partieAvion(k).Icm * transpose(R);
       
    end 

    ux = R*ux;
end

names = {partieAvion.name};                                         
centreMasseMoteurDroit = partieAvion(strcmp(names,"moteur_droit")).r;
centreMasseMoteurGauche = partieAvion(strcmp(names,"moteur_gauche")).r;
centreMasseAileDroite = partieAvion(strcmp(names,"aile_droite")).r;
centreMasseAileGauche = partieAvion(strcmp(names,"aile_gauche")).r;

   centreMasseMoteurDroitApplique = centreMasseMoteurDroit - (LONGUEUR_MOTEUR/2)*ux;
   centreMasseMoteurGaucheApplique = centreMasseMoteurGauche - (LONGUEUR_MOTEUR/2)*ux;
   centreMasseForcePorteeApplique = 0.5 *(centreMasseAileDroite + centreMasseAileGauche);
   centreMasseForcePorteeApplique = [centreMasseForcePorteeApplique(1); centreMasseForcePorteeApplique(2);0 ];

   pts = struct('pos',{centreMasseMoteurDroitApplique, centreMasseMoteurGaucheApplique, centreMasseForcePorteeApplique});  


   forceMoteurDroit = Forces(1)*ux;
   forceMoteurGauche = Forces(2)*ux;
   ForcePorte  = [0;0;Forces(3)];
   F =  [ forceMoteurDroit, forceMoteurGauche, ForcePorte];



   MI = calculInertie(partieAvion,pcm);
   momentTotal = calculMomentForce(pts, F, pcm);
   momentCinetique = MI * va;
   aa = MI\( momentTotal + cross(momentCinetique, va));

end


 

function [momentTotal] = calculMomentForce(vecteurPosForce,F, pcm)
    momentTotal = zeros(3,1);
    for k = 1:numel(vecteurPosForce)
        momentTotal = momentTotal + cross( (vecteurPosForce(k).pos - pcm) , F(:,k));
    end

end 


function [pcm, partieAvion] = calculCentreMasse(posA)

baseCabine =  [-(3/4)*HAUTEUR_CABINE ; 0 ; 0 ];
cabineCentreMasse = posA + baseCabine;


fuselageCentreMasse = posA + [-(HAUTEUR_CABINE + (LONGUEUR_FUSELAGE/2));0 ;0];

ailesCentreMasseDroite = [CENTRE_X_AILE ; CENTRE_Y_AILE ; EPAISSEUR_AILE];
ailesCentreMasseGauche = [CENTRE_X_AILE ; -CENTRE_Y_AILE ; EPAISSEUR_AILE];

moteurCentreMasseDroite = [5 ;RAYON_FUSELAGE+RAYON_MOTEUR; RAYON_FUSELAGE + EPAISSEUR_AILE];
moteurCentreMasseGauche = [5 ; -(RAYON_FUSELAGE+RAYON_MOTEUR); RAYON_FUSELAGE + EPAISSEUR_AILE];

AileronCentreMasse = [(LARGEUR_AILERON/2); 0 ; ((2*RAYON_FUSELAGE)+ EPAISSEUR_AILERON) ] ;

Sum = (cabineCentreMasse * MASSE_CABINE) + (ailesCentreMasseDroite * MASSE_AILE)  ...
    + (ailesCentreMasseGauche * MASSE_AILE) + (moteurCentreMasseDroite * MASSE_MOTEUR )...  
    + (moteurCentreMasseGauche  * MASSE_MOTEUR) +(AileronCentreMasse * MASSE_AILERON);

masseTotale = (MASSE_CABINE) + 2*(MASSE_AILE) + 2*(MASSE_MOTEUR) + (MASSE_AILERON);

pcm = Sum/masseTotale;


partieAvion = [ ...
    struct('name',"cabine",        'm',MASSE_CABINE,   'r',cabineCentreMasse,        'Icm',INERTIE_CENTRE_CABINE), ...
    struct('name',"fuselage",      'm',MASSE_FUSELAGE, 'r',fuselageCentreMasse,      'Icm',INTERTIE_CENTRE_FUSELAGE), ...
    struct('name',"aile_droite",   'm',MASSE_AILE,     'r',ailesCentreMasseDroite,   'Icm',INTERTIE_CENTRE_AILE), ...
    struct('name',"aile_gauche",   'm',MASSE_AILE,     'r',ailesCentreMasseGauche,   'Icm',INTERTIE_CENTRE_AILE), ...
    struct('name',"moteur_droit",  'm',MASSE_MOTEUR,   'r',moteurCentreMasseDroite,  'Icm',INTERTIE_CENTRE_MOTEUR), ...
    struct('name',"moteur_gauche", 'm',MASSE_MOTEUR,   'r',moteurCentreMasseGauche,  'Icm',INTERTIE_CENTRE_MOTEUR), ...
    struct('name',"aileron",       'm',MASSE_AILERON,  'r',AileronCentreMasse,       'Icm',INTERTIE_CENTRE_AILERON) ...
];



end 


function [inertieTotal] = calculInertie(partieAvion, pcm)
    inertieTotal = zeros(3,3);  
    for k = 1:numel(partieAvion)
        r = partieAvion(k).r;
        d =  pcm - r;
        dx = d(1); 
        dy = d(2); 
        dz = d(3);
        T =[ dy^2 + dz^2,  -dx*dy,      -dx*dz;
              -dx*dy,       dx^2 + dz^2, -dy*dz;
              -dx*dz,       -dy*dz,      dx^2 + dy^2 ];
    
        inertieTotal = inertieTotal + partieAvion(k).Icm + partieAvion(k).m * T;
    end
end
