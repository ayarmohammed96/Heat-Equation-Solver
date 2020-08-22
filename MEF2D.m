function [Kb,Fb,U] = MEF2D(f,elem,k,n,X,T,X1,T1,b,Ngauss)
% fonction traite l'equation alpha*laplace(u)=f sur Omega avec u|Omega=0
% moyennant la méthode des éléments finis de type Galerkin Pk
%Entrées
%          f:  fonction second membre
%       elem:  type d'élément (0 pour quadrilatères et 1 pour triangles)
%          k:  type de la méthodes des éléments finis (Pk)
%          X:  table des coordonnées
%          T:  table de connectivité
%          b:  noeuds sur le bords
%     Ngauss:  Nombre des noeuds de Gauss dans l'élément courant 
%Sorties:
%          U:  Vecteurs solutions aux noeuds du maillage.
%
% --------Quadrature de Gauss-----------
[Xgauss,Wgauss]=Quadrature(elem, Ngauss);
%
%---------Fonctions chapeaux------------
[phi,grad] = FoncChap(elem,k,n,Xgauss);
%
%-------- Assemblage--------------------
[ K,F ] = Assemblage2D(elem,k,n,X,T,X1,T1,f,phi,grad,Xgauss,Wgauss,Ngauss);
%
%--------Conditions aux bords-----------
%On met "0" dans les lignes et les colonnes qui correspondent aux noeuds
%bords dans K et F
K(b,:)=0; K(:,b)=0; 
F(b)=0; 
%On met "1" sur la diagonale de K correspodant aux aux noeuds bords 
K(b,b)=speye(length(b),length(b)); % cas où K est sparse
%K(b,b)=eye(length(b),length(b));
% K et F après conditions aux bords
Kb=K; Fb=F; 
%
U=Kb\Fb;
end

