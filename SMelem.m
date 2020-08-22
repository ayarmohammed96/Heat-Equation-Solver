function [ Fe ] = SMelem(f,elem,k,n,X,T,ie,Ngauss)
%calcul du second membre élémentaire Fe
%
[Xgauss,Wgauss]=Quadrature(elem, Ngauss);
%fonctions chapeaux
phi = FoncChap(elem,k,n,Xgauss);
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
%

if k==1
Fe=zeros(3*k,1); %initialiser les tailles de Me et Fe 
%      Noeuds du triangle courant (sommet dans le cas où Ngeom=3)
x=X(T(ie,1:Ngeom),:); 
%      matice du changement de variable
J=[(x(2,:)-x(1,:))',(x(3,:)-x(1,:))'];

%
%On calcule Xc=X1+J*S aux noeuds de Gauss dont le nombre est Ngauss
  Xc=zeros(2,Ngauss);
  for i=1:Ngauss      
        Xc(:,i)=x(1,:)'+(J*Xgauss(i,:)');
  end  
%
%on calcule Fe1 Fe2 et Fe2 les composantes de Fe
  for i=1:3
      s=0;
      for j=1:Ngauss
       s=s+ (Wgauss(j)*f(Xc(1,j),Xc(2,j))*phi(j,i));
      end
      Fe(i)=det(J)*s;
  end

else 
        error(' ce type Pk n est pas encore programmé') 
end
end

