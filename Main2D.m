clear all, close all, home
%------------------------
% Les entrées
%------------------------
elem=input('Choisir elment(0 quadrilatère, 1 triangle) :  ');
n=input('Choisir 1 pour maillage non conforme 0 pour conforme :  ' );
while (n~=0)&&(n~=1)
clear n;
n=input('nombre choisi ne correspond a aucun type de maillage, saisissez 1 pour non conforme ou 0 pour conforme :  ' );
end
k=input('Choisir k pour type Pk :  ');
Ngauss=input('Choisir le nombre de points de Gauss :  ');
npx=input('Entrer nombre de points suivant (Ox) :  ');
npy=input('Entrer nombre de points suivant (Oy) :  ');
%------------------------------
% Les données du maillage
x1=-2;x2=2;
y1=-1;y2=1;
%
% Second membre
%
f=@(x,y) -2*x^2-2*y^2+10;
% Solution exacte
ue=@(x,y) (x^2-4)*(y^2-1);
%--------------------------------
% maillage
[X,T,X1,T1,b] = Maillage2d(elem,k,n,x1,x2,y1,y2,npx,npy);
% calcul de la solution
[Kb,Fb,U] = MEF2D(f,elem,k,n,X,T,X1,T1,b,Ngauss);
%------------------------------
% Affichage et comparaison
if n==0
% Solution approchée   
figure(2); clf;
trisurf(T,X(:,1),X(:,2),0*X(:,1),U,'edgecolor','k','facecolor','interp');
view(2),axis([x1 x2 y1 y2]),axis equal,colorbar;
figure(3); clf;
trisurf(T,X(:,1),X(:,2),U);
%}
% Solution exacte
figure(4); clf;
uee=(X(:,1).^2-4).*(X(:,2).^2-1);
trisurf(T,X(:,1),X(:,2),0*X(:,1),uee,'edgecolor','k','facecolor','interp');
view(3),axis([x1 x2 y1 y2]),axis equal,colorbar;
figure(5); clf;
trisurf(T,X(:,1),X(:,2),uee);
elseif n==1
% Solution approchée
figure(2); clf;
trisurf(T1,X1(:,1),X1(:,2),0*X1(:,1),U,'edgecolor','k','facecolor','interp');
view(2),axis([x1 x2 y1 y2]),axis equal,colorbar;
figure(3); clf;
trisurf(T1,X1(:,1),X1(:,2),U);
%}
% Solution exacte
figure(4); clf;
uee=(X1(:,1).^2-4).*(X1(:,2).^2-1);
trisurf(T1,X1(:,1),X1(:,2),0*X1(:,1),uee,'edgecolor','k','facecolor','interp');
view(3),axis([x1 x2 y1 y2]),axis equal,colorbar;
figure(5); clf;
trisurf(T1,X1(:,1),X1(:,2),uee);
%--------------------------------
end
% Calcul de l'erreur en norme L2
error_L2=norm((U-uee),2);
err_rela=error_L2/norm(uee,2);
fprintf('\t\t %4d \t\t %20.16e \t\t %20.16e\n\n', npx*npy,error_L2,err_rela);
fprintf('On a fini normalement\n');
