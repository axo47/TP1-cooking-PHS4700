% Script de test pour Devoir1
% Cas 1
posCas1=[3.82+22.95;0;1.595];
arCas1=0.0;
vaCas1=[0;0;0];
ForcesCas1=[11e6; 11e6; 260e6];
[pcmCas1 MICas1 aaCas1]=Devoir1(posCas1,arCas1,vaCas1,ForcesCas1);
fprintf('\nCas 1\n');
fprintf('  Position avion (m) = ( %10.5f,  %10.5f,  %10.5f )\n',posCas1(1),posCas1(2),posCas1(3));
fprintf('  Rotation avion (r) = %10.5f\n',arCas1);
fprintf('  Vitesse angulaire avion (r/s) = ( %10.5f,  %10.5f,  %10.5f )\n',vaCas1(1),vaCas1(2),vaCas1(3));
fprintf('  Forces (N) = ( %10.5f,  %10.5f,  %10.5f ) \n',ForcesCas1(1),ForcesCas1(2),ForcesCas1(3));
fprintf('\nResultats avion sur la piste \n');
fprintf('  Centre de masse (m) = ( %10.5f,  %10.5f,  %10.5f )\n',pcmCas1(1),pcmCas1(2),pcmCas1(3));
fprintf('  Moment inertie (kg/m^2) =\n   %14.0f & %14.0f & %14.0f  \n   %14.0f & %14.0f & %14.0f  \n   %14.0f & %14.0f & %14.0f \n',...
    MICas1(1,1),MICas1(1,2),MICas1(1,3),MICas1(2,1),MICas1(2,2),MICas1(2,3),MICas1(3,1),MICas1(3,2),MICas1(3,3));
fprintf('  acc angulaire (r/s^2) = ( %10.5f,  %10.5f,  %10.5f ) \n\n',aaCas1(1),aaCas1(2),aaCas1(3));

% Cas 2
posCas2=[4198;0;618];
arCas2=-0.15;
vaCas2=[0.0; -0.003; -0.01];
ForcesCas2=[8e6; 11e6; 260e6];
[pcmCas2 MICas2 aaCas2]=Devoir1(posCas2,arCas2,vaCas2,ForcesCas2);
fprintf('\nCas 2\n');
fprintf('  Position avion (m) = ( %10.5f,  %10.5f,  %10.5f )\n',posCas2(1),posCas2(2),posCas2(3));
fprintf('  Rotation avion (r) = %10.5f \n',arCas2);
fprintf('  Vitesse angulaire avion (r/s)= ( %10.5f,  %10.5f,  %10.5f )\n',vaCas2(1),vaCas2(2),vaCas2(3));
fprintf('  Forces (N) = ( %10.5f,  %10.5f,  %10.5f )\n',ForcesCas2(1),ForcesCas2(2),ForcesCas2(3));
fprintf('\nResultats avion en vol\n');
fprintf('  Centre de masse (m) = ( %10.5f,  %10.5f,  %10.5f )\n',pcmCas2(1),pcmCas2(2),pcmCas2(3));
fprintf('  Moment inertie  (kg/m^2) =\n   %14.0f & %14.0f & %14.0f \\\\ \n   %14.0f & %14.0f & %14.0f \\\\ \n   %14.0f & %14.0f & %14.0f \\\\ \n',...
    MICas2(1,1),MICas2(1,2),MICas2(1,3),MICas2(2,1),MICas2(2,2),MICas2(2,3),MICas2(3,1),MICas2(3,2),MICas2(3,3));
fprintf('  acc angulaire (r/s^2) = ( %10.5f,  %10.5f,  %10.5f )\n\n',aaCas2(1),aaCas2(2),aaCas2(3));