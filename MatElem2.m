function [Ke] = MatElem2(elem,k,X,T,X1,T1,ie,phi,grad,Xgauss,Wgauss)
%Fonction qui calcule les matrices élémentaires 
%Entrées:
%          k: type de l'élément fini (Pk)
%          X: table des coordonnées
%          T: tbale des connectivités
%          ie: indice de l'éléménts courant
%          f: fonction du second membre
%         phi: fonctions chapeaux dans TR
%        grad: gradiants des fonctions chapeaux dans TR
%      Xgauss: les noeuds de Gauss dans l'élément de Référence
%      Wgauss: les poids de Gauss dans l'élément de Référence
%Sorties:
%          Ke: matrice élémentaire (des gradiants)
%          Fe: second membre élémentaire
% 
if elem == 0
    if k == 1
        Ngeom = 4;
    elseif k == 2
        Ngeom=9;
    end
elseif elem == 1  
    if k == 1
        Ngeom = 3;             
    elseif k == 2 
        Ngeom = 6;
    elseif k == 3
        Ngeom = 9;
    end
end
if k==1 
%      Noeuds du triangle courant (sommet dans le cas où Ngeom=3)
x=X(T(ie,1:Ngeom),:);

%      matice du changement de variable
J=[(x(2,:)-x(1,:))',(x(3,:)-x(1,:))'];

%       surface de l'élément de référence
     air=0.5*abs(det([0 1 0;0 0 1;1 1 1]));
%  
   KKe=zeros(3*k);%initialiser les tailles de Ke
%
%calcul de la matrice élémentaire Ke
   for i=1:3 
     for j=i:3
        KKe(i,j)=(J\grad(:,i))'*(J\grad(:,j));
        KKe(j,i)=KKe(i,j);
     end
   end
   Ke=air*det(J)*KKe;   
 %
else 
        error(' ce type Pk n est pas encore programmé') 
end  
end

