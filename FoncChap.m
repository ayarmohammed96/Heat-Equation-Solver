function [phi,grad] = FoncChap(elem,k,n,Xgauss)
%cette fonction donne les fonctions de formes dans l'element de référence
%        N3=(0,1)
%               | \
%               |  \
%               |   \
%               | TR \
%               |_ _ _\
%        N1=(0,0)      (1,0)=N2
%Entrées:
%            k: type des éléments finis Pk
%       Xgauss: Noeuds de Gauss dans TR
%Sorties:
%         phi: tableau dont les colonnes contiennent les fonctions chapeaux 
%              phi1, phi2 et phi3 dans cet ordre
%        grad: tableau dont les colonnes contiennent les gradiants  
%              des fonctions chapeaux phi1, phi2 et phi3 dans cet ordre
%
s = Xgauss(:,1); 
t = Xgauss(:,2); 
if  elem==1
      if n==0
        if k == 1  %P1
        phi   = [1-(s+t),s,t]; 
        grad  = [-ones(2,1),eye(2)];
        end
      elseif n==1
           if k==1
           phi   = [-2*t+1,-1+2*(s+t),-2*s+1]; 
           grad  = [0 2 -2;-2 2 0];
           elseif k==2   
           phi   = [-2*t+1,-1+2*(s+t),-2*s+1]; 
           grad  = [0 2 -2;-2 2 0];
           end  
       end
elseif elem==0
error('non programmé');
end
end

