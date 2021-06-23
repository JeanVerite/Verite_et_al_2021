clc 
clear all

% ------------------------------------------------------------------------------------------------------------------------------------------

%        IMPORTATION DES DONNEES .TXT DE DES FICHIERS VECTEURS

% Ouverture fichier .txt pour triangulation
[filename2, pathname2, filterindex2] = uigetfile('*.txt', 'fichier vecteurs vitesse'); % Chemin d'accès du premier fichier
nb = dir (pathname2);
long = length (find([nb.isdir]==0));
premier_fichier=fopen([num2str(pathname2),num2str(filename2)]); % Lecture du fichier
a=fscanf(premier_fichier,'%f',[4,inf]);
X1=a(1,:);
Y1=a(2,:);
V_x=a(3,:);
V_y=a(4,:);
% pas_temps=input('Pas de temps en secondes');
pas_temps=90;
[aa,bb]=size(X1);
X2=ones(1,bb);
Y2=ones(1,bb);
for i=1:bb;
    X2(i,1)=X1(i)+V_x(i)*pas_temps;
    Y2(i,1)=Y1(i)+V_y(i)*pas_temps;
end;

% ------------------------------------------------------------------------------------------------------------------------------------------------

%                      TRIANGULATION DE DELAUNAY 

DT1 = delaunay(X1,Y1); % triangulation pour positions à t2
sz2=size(DT1,1); % nb de triangle créés

% ------------------------------------------------------------------------------------------------------------------------------------------------

%       CREATION DES FICHIERS TXT D'ENREGISTREMENT DES DONNEES

numero_fichier1=input('Nom du fichier : ');
fgv1=([num2str(pathname2),'Données Déformation',num2str(numero_fichier1),'.txt']); % Fichiers .txt 
fid1=fopen(fgv1,'wt');

% ---------------------------------------------------------------------------------------------------------------

% BOUCLE CALCUL : ROTATION ELLISPE DE DEFORMATION + VITESSE DE DEFORMATION

for i=1:sz2;  % POUR CHAQUE TRIANGLE
    
    % COORDONNEES DU TRIANGLE ABC
    
        % Calcul des coordonnées x et y des vecteurs A et B à t1 et t2
    
    vecteur1(i,1)=X1(DT1(i,2))-X1(DT1(i,1)); % xA à t1
    vecteur2(i,1)=X2(DT1(i,2))-X2(DT1(i,1)); % xA à t2
    vecteur1(i,2)=Y1(DT1(i,2))-Y1(DT1(i,1)); % yA à t1
    vecteur2(i,2)=Y2(DT1(i,2))-Y2(DT1(i,1)); % yA à t2
    vecteur1(i,3)=X1(DT1(i,3))-X1(DT1(i,1)); % xB à t1
    vecteur2(i,3)=X2(DT1(i,3))-X2(DT1(i,1)); % xB à t2
    vecteur1(i,4)=Y1(DT1(i,3))-Y1(DT1(i,1)); % yB à t1
    vecteur2(i,4)=Y2(DT1(i,3))-Y2(DT1(i,1)); % yB à t2
    
        % Calcul coordonnées centre de gravité du triangle à t2
    
    centre_triangle_x=(X2(DT1(i,1))+X2(DT1(i,2))+X2(DT1(i,3)))/3; % coordonnées en x du centre de gravité du triangle
    centre_triangle_y=(Y2(DT1(i,1))+Y2(DT1(i,2))+Y2(DT1(i,3)))/3; % coordonnées en y du centre de gravité du triangle
    
    
    % COMPOSANTES VITESSE Vx et Vy entre t2  et t1 pour A, B, C
    
        % Point 1 (A)
    vitesse_x_1=(X1(DT1(i,1))-X2(DT1(i,1)))/pas_temps;
    vitesse_y_1=(Y1(DT1(i,1))-Y2(DT1(i,1)))/pas_temps;
        % Point 2 (B)
    vitesse_x_2=(X1(DT1(i,2))-X2(DT1(i,2)))/pas_temps;
    vitesse_y_2=(Y1(DT1(i,2))-Y2(DT1(i,2)))/pas_temps;
        % Point 3 (C)
    vitesse_x_3=(X1(DT1(i,3))-X2(DT1(i,3)))/pas_temps;
    vitesse_y_3=(Y1(DT1(i,3))-Y2(DT1(i,3)))/pas_temps;
    
    % Calcul aire triangle
    
    aire_triangle=0.5*(X2(DT1(i,2)).*Y2(DT1(i,3))-X2(DT1(i,3)).*Y2(DT1(i,2))+Y2(DT1(i,1)).*X2(DT1(i,3))-X2(DT1(i,1)).*Y2(DT1(i,3))+X2(DT1(i,1)).*Y2(DT1(i,2))-X2(DT1(i,2)).*Y2(DT1(i,1)));
    
    % TENSEUR DE DEFORMATION
    
    dvx_x=(1/(2*aire_triangle))*(vitesse_x_1*(Y2(DT1(i,2))-Y2(DT1(i,3)))+vitesse_x_2*(Y2(DT1(i,3))-Y2(DT1(i,1)))+vitesse_x_3*(Y2(DT1(i,1))-Y2(DT1(i,2))));
    dvx_y=(1/(2*aire_triangle))*(vitesse_x_1*(X2(DT1(i,3))-X2(DT1(i,2)))+vitesse_x_2*(X2(DT1(i,1))-X2(DT1(i,3)))+vitesse_x_3*(X2(DT1(i,2))-X2(DT1(i,1))));
    dvy_x=(1/(2*aire_triangle))*(vitesse_y_1*(Y2(DT1(i,2))-Y2(DT1(i,3)))+vitesse_y_2*(Y2(DT1(i,3))-Y2(DT1(i,1)))+vitesse_y_3*(Y2(DT1(i,1))-Y2(DT1(i,2))));
    dvy_y=(1/(2*aire_triangle))*(vitesse_y_1*(X2(DT1(i,3))-X2(DT1(i,2)))+vitesse_y_2*(X2(DT1(i,1))-X2(DT1(i,3)))+vitesse_y_3*(X2(DT1(i,2))-X2(DT1(i,1))));
    
    Exx=dvx_x;
    Exy=0.5*(dvx_y+dvy_x);
    Eyy=dvy_y;
    
    % ELLIPSE DE DEFORMATION
    
        % Création matrices carrées des sets de vecteurs A-B à t1 et t2
    
    A2=[vecteur2(i,1),vecteur2(i,3);vecteur2(i,2),vecteur2(i,4)];
    A1=[vecteur1(i,1),vecteur1(i,3);vecteur1(i,2),vecteur1(i,4)];
    INV1=inv(A1);
    
        % Calcul matrice de déformation D entre t2-t1 
    
    D=A2*INV1; % D (a b ; c d)
    vecteurD(i,1)=D(1,1); % a
    vecteurD(i,2)=D(1,2); % b
    vecteurD(i,3)=D(2,1); % c
    vecteurD(i,4)=D(2,2); % d
   
        % Résolution équation du 2nd degré pour retrouver m
    
            % u = (1 ; m)  v = (-m ; 1) | u' = D*u  v' = D*v    u'T.v'=0
            % D(a b ; c d) | m²(-ab-cd)+m(-a²+b²-c²+d²)+ab+cd = 0
   
    grandA=(-vecteurD(i,1)*vecteurD(i,2))-(vecteurD(i,3)*vecteurD(i,4)); % A = -ab-cd
    grandB=-(vecteurD(i,1)^2)+ vecteurD(i,2)^2-(vecteurD(i,3)^2)+vecteurD(i,4)^2; % B = -a²+b²-c²+d²
    grandC=vecteurD(i,1)* vecteurD(i,2)+vecteurD(i,3)*vecteurD(i,4); % C = ab+cd
    
    delta=grandB^2-(4*grandA*grandC);
    m(i,1)=(-grandB+sqrt(delta))/(2*grandA);
    m(i,2)=(-grandB-sqrt(delta))/(2*grandA);
    valeurM=max(m(i,:));
    
        % Définition coordonnées des vecteurs U-V et U'-V' 
    
    x_vecteurU=1;                   %Ux
    y_vecteurU=valeurM;             %Uy
    x_vecteurV=-(valeurM);          %Vx
    y_vecteurV=1;                   %Vy
    
    x_vecteurU_prim=vecteurD(i,1)+vecteurD(i,2)*valeurM;      %U'x
    y_vecteurU_prim=vecteurD(i,3)+vecteurD(i,4)*valeurM;      %U'y
    x_vecteurV_prim=-(vecteurD(i,1)*valeurM)+vecteurD(i,2);   %V'x
    y_vecteurV_prim=-(vecteurD(i,3)*valeurM)+vecteurD(i,4);   %V'y
    
        % Calcul de l'angle de rotation des vecteurs U et V de l'ellipse de
        % déformation (+ NoData si angle aberrant)
    
    cos_teta=(x_vecteurU*x_vecteurU_prim+y_vecteurU*y_vecteurU_prim)/(sqrt(x_vecteurU^2+y_vecteurU^2)*sqrt(x_vecteurU_prim^2+y_vecteurU_prim^2));
    if acos(cos_teta)*180/pi > 20;
        teta = 0;
    else teta = acos(cos_teta)*180/pi;
    end ;
    
    TF=isnan(teta);
    
%         % Calcul longueurs demi-axes de l'ellipse et rapport U'/U et V'/V
    longU=sqrt((x_vecteurU)^2+(y_vecteurU)^2);
    longUprim=sqrt((x_vecteurU_prim)^2+(y_vecteurU_prim)^2);
    longV=sqrt((x_vecteurV)^2+(y_vecteurV)^2);
    longVprim=sqrt((x_vecteurV_prim)^2+(y_vecteurV_prim)^2);
   
    r_U_Uprim=(longUprim/longU)-1;
    if r_U_Uprim>0;
        elongU=1;
    else elongU=0;
    end;
    r_V_Vprim=(longVprim/longV)-1;
    if r_V_Vprim>0;
        elongV=1;
    else elongV=0;
    end;
%     
%     % Calcul coordonnées vecteurs "différentiel d'allongement de
%     % l'ellipse"
%     
     x_bis_Uprim=x_vecteurU_prim.*r_U_Uprim/longUprim;

     y_bis_Uprim=y_vecteurU_prim.*r_U_Uprim/longUprim;
   
     x_bis_Vprim=x_vecteurV_prim.*r_V_Vprim/longVprim;

     y_bis_Vprim=y_vecteurV_prim.*r_V_Vprim/longVprim;
    
    % Remplissage matrice data finales
    
        % Position centre de gravité triangle
    mat_def(i,1)=centre_triangle_x;
    mat_def(i,2)=centre_triangle_y;
    mat_def(i,3)=teta;
    mat_def(i,4)=Exy;
    
%         % Rapport vecteur U'/U 
    mat_def(i,3)=x_bis_Uprim.*400; % *400 pour dilater la longueur du vecteur sur le graphique
    mat_def(i,4)=y_bis_Uprim.*400;
%     
%         % Rapport vecteur V'/V 
    mat_def(i,5)=x_bis_Vprim.*400;
    mat_def(i,6)=y_bis_Vprim.*400;
% 
%         % Rotation axes de l'ellispoide entre t1 et t2
    mat_def(i,7)=teta;
    mat_def(i,8)=elongU;
    mat_def(i,9)=elongV;
%     
%         % Tenseur déformation
    mat_def(i,10)=Exx;
    mat_def(i,11)=Exy;
    mat_def(i,12)=Eyy;
 
  if mat_def(i,3) ~= 0;
    if mat_def(i,4) ~= 0;
        if TF == 0;
            fprintf(fid1,'%f %f %f %f\n',mat_def(i,1),mat_def(i,2),mat_def(i,3),mat_def(i,4));
        end;
    end;
  end;
   
U_prim_demi_grand_axe(i)=sqrt((mat_def(i,3))^2+(mat_def(i,4))^2);
V_prim_demi_grand_axe(i)=sqrt((mat_def(i,5))^2+(mat_def(i,6))^2);
       
end;

% ----------------------------------------------------------------------------------------------------------------

% % %                           SORTIES GRAPHIQUES
% 
% % Triangulation de Delaunay
f2=figure;
figure(f2);
% triplot(DT1,X2,Y2); 
hold on
% 
% % Allongement du demi-grand axe pour ellispoide de déformation à t2
% 
for i=1:sz2; % pour chaque triangle
    if U_prim_demi_grand_axe(i)>V_prim_demi_grand_axe(i); % si U'>V'
        if mat_def(i,8)==1; % et que l'allongement de U' est positif
        quiver(mat_def(i,1),mat_def(i,2),mat_def(i,3),mat_def(i,4),'k'); % tracer le vecteur allongement de U'
        end;
    else; % si U'<V'
        if mat_def(i,9)==1; % et que l'allongement de V' est positif
        quiver(mat_def(i,1),mat_def(i,2),mat_def(i,5),mat_def(i,6),'b'); % tracer le vecteur allongement de U'
        end;
    end;
end;

hold off

% ----------------------------------------------------------------------------------------------------------
