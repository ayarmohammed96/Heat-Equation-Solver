function [ K,F ] = Assemblage2D(elem,k,n,X,T,X1,T1,f,phi,grad,Xgauss,Wgauss,Ngauss)
% Assemblage des matrices �l�mentaires "Ke" dans la matrice globale K
% Assemblage des seconds membres �l�mentaires "Fe" dans le second... 
% membre global F
% On utilise la fonction "MatElem(k,X,T,ie,f,phi,grad,Xgauss,Wgauss)"
% 
%Entr�es:
%       elem:  type d'�l�ment (0 pour quadrilat�res et 1 pour triangles)
%          X: table des coordonn�es
%          T: table des connectivit�s
% 

if n==0
Nn=size(X,1);  % nombre des noeuds
Nt=size(T,1);  % nombre des �l�ments
K=sparse(zeros(Nn,Nn));% initialisation de K comme matrice creuse
%K=zeros(Nn,Nn);% initialisation de K 
F=zeros(Nn,1);% initialisation de F comme matrice creuse
for ie = 1:Nt 
    Tie = T(ie,:);
    [Ke] = MatElem2(elem,k,X,T,X1,T1,ie,phi,grad,Xgauss,Wgauss);
    [Fe] = SMelem(f,elem,k,n,X,T,ie,Ngauss);
    K(Tie,Tie)=K(Tie,Tie)+Ke;
    F(Tie)=F(Tie)+Fe;
    clear Ke Fe; 
end
else
Nn=size(X1,1);  % nombre des noeuds
Nt=size(T1,1);  % nombre des �l�ments
K=sparse(zeros(Nn,Nn));% initialisation de K comme matrice creuse
%K=zeros(Nn,Nn);% initialisation de K 
F=zeros(Nn,1);% initialisation de F comme matrice creuse   
for ie = 1:Nt 
    Tie = T1(ie,:);
    [Ke] = MatElem2(elem,k,X,T,X1,T1,ie,phi,grad,Xgauss,Wgauss);
    [Fe] = SMelem(f,elem,k,n,X,T,ie,Ngauss);
    K(Tie,Tie)=K(Tie,Tie)+Ke;
    F(Tie)=F(Tie)+Fe;
    clear Ke Fe; 
end
end
K=sparse(K);
end

