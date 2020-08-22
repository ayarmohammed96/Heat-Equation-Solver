function [Ke] = MatElem2(elem,k,X,T,X1,T1,ie,phi,grad,Xgauss,Wgauss)
%Fonction qui calcule les matrices �l�mentaires 
%Entr�es:
%          k: type de l'�l�ment fini (Pk)
%          X: table des coordonn�es
%          T: tbale des connectivit�s
%          ie: indice de l'�l�m�nts courant
%          f: fonction du second membre
%         phi: fonctions chapeaux dans TR
%        grad: gradiants des fonctions chapeaux dans TR
%      Xgauss: les noeuds de Gauss dans l'�l�ment de R�f�rence
%      Wgauss: les poids de Gauss dans l'�l�ment de R�f�rence
%Sorties:
%          Ke: matrice �l�mentaire (des gradiants)
%          Fe: second membre �l�mentaire
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
%      Noeuds du triangle courant (sommet dans le cas o� Ngeom=3)
x=X(T(ie,1:Ngeom),:);

%      matice du changement de variable
J=[(x(2,:)-x(1,:))',(x(3,:)-x(1,:))'];

%       surface de l'�l�ment de r�f�rence
     air=0.5*abs(det([0 1 0;0 0 1;1 1 1]));
%  
   KKe=zeros(3*k);%initialiser les tailles de Ke
%
%calcul de la matrice �l�mentaire Ke
   for i=1:3 
     for j=i:3
        KKe(i,j)=(J\grad(:,i))'*(J\grad(:,j));
        KKe(j,i)=KKe(i,j);
     end
   end
   Ke=air*det(J)*KKe;   
 %
else 
        error(' ce type Pk n est pas encore programm�') 
end  
end

