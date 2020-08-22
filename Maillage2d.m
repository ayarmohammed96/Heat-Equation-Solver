function [X,T,X1,T1,b] = Maillage2d(elem,k,n,x1,x2,y1,y2,npx,npy)
% crée un maillage uniforme sur un domaine 
% rectangulaire [x1,x2]x[y1,y2]
%
% Entrées:    
%   x1,x2,y1,y2:    coordonnées des coins 
%   npx,npy:        nombres de noeuds dans chaque direction 
% Sorties:   
%          X:  table des coordonnées
%          T:  table de connectivité
%          b:  noeuds sur le bords

% Nombre des éléments dans chaque direction
 nx=npx-1 ;  ny=npy-1 ;
 X = zeros((npx)*(npy),2);
 
%
hx = (x2-x1)/nx;
hy = (y2-y1)/ny;
xs = linspace(x1,x2,npx)'; 
unos = ones(npx,1);
%
% Coordonnées des noeuds
%P1
yys = linspace(y1,y2,npy);
for i=1:npy
    ys = yys(i)*unos; 
    posi = (i-1)*(npx)+1:i*(npx); 
    X(posi,:)=[xs,ys];
end
% Coordonnées des noeuds 
%NP1
if n==1
    if k==1
    %NP1
X1 = zeros((2*npx-1)*(npy-1)+(nx*npy),2);
yys1 = linspace(y1,y2,npy);
yys2 = linspace(y1+hy/2,y2-hy/2,ny);
indice1=0;indice2=0;
for i=1:2*npy-1
    if mod(i,2)==1
    xs1 = linspace(x1+hx/2,x2-hx/2,nx)';
    ys1 = yys1(i/2+1/2)*ones(nx,1);
    posi = indice2+1:indice2+nx;
    X1(posi,:)=[xs1,ys1];
    indice1= max(posi);
    else
    xs2 = linspace(x1,x2,2*npx-1)';
    ys2 = yys2(i/2)*ones(2*npx-1,1); 
    posi = indice1+1:indice1+(2*npx-1);
    X1(posi,:)=[xs2,ys2];
    indice2= max(posi);
    end    
end
    elseif k==2
    %NP2
    X1 = zeros(2*(2*npx-1)*ny+(2*nx*npy),2);
yys1 = linspace(y1,y2,npy);
yys2 = y1:hy/3:y2;
yys2([1:3:end])=[];
indice1=0;indice2=0;indice3=0;
for i=1:npy+2*ny
   xs1 = linspace(x1,x2,npx+2*nx)';
   xs1([1:3:npx+2*nx])=[];
   xs21 = linspace(x1,x2,npx+2*nx)';
   xs21([2:3:npx+2*nx])=[];
   xs22 = linspace(x1,x2,npx+2*nx)';
   xs22([3:3:npx+2*nx])=[];
   if mod(i-1,3)==0
        if mod(i,2)==1
    ys1 = yys1(i/3+2/3)*ones(2*nx,1);
    posi = indice3+1:indice3+2*nx;
    X1(posi,:)=[xs1,ys1];
    indice1= max(posi);
        else
    ys1 = yys1(i/3+2/3)*ones(2*nx,1);
    posi = indice3+1:indice3+2*nx;
    X1(posi,:)=[xs1,ys1];
    indice1= max(posi);    
        end
   else
       if mod(i+1,3)==0
    ys2 = yys2((2*i)/3-1/3)*ones(2*npx-1,1); 
    posi = indice1+1:indice1+(2*npx-1);
    X1(posi,:)=[xs21,ys2];
    indice2= max(posi);
       else
    ys2 = yys2((2*i)/3)*ones(2*npx-1,1); 
    posi = indice2+1:indice2+(2*npx-1);
    X1(posi,:)=[xs22,ys2];
    indice3= max(posi);    
       end
    end    
end
    end
end
% P1
% Connectivité
if elem==0
T = zeros(nx*ny,4);    % Maillage Quadrilatère
        for i=1:ny
            for j=1:nx
                ielem = (i-1)*nx+j;
                inode = (i-1)*(npx)+j;
                T(ielem,:) = [inode   inode+1   inode+npx inode+npx+1];
               
            end   
        end
elseif elem==1         % Maillage triangulaire
    
    T = zeros(nx*ny,3);
        for i=1:ny
            for j=1:nx
                ielem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*(npx)+j;
                T(ielem,:) = [inode   inode+1   inode+(npx)];
                T(ielem+1,:) = [inode+1+npx   inode+npx inode+1];
                
            end   
        end
     %}   
    if n==1
        if k==1
        %NP1
    T1 = zeros(nx*ny,3);
        for i=1:ny
            for j=1:nx
                ielem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*(2*npx-1+nx)+j;
                T1(ielem,:) = [inode inode+(nx+j) inode+nx+j-1];
                T1(ielem+1,:) = [inode+2*npx-1+nx inode+(nx+j) inode+(nx+j)+1  ] ;
             
            end
        end
        elseif k==2
        %NP2
     T1 = zeros(nx*ny,3);
        for i=1:ny
            for j=1:nx
                ielem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*(npx)+j;
                T(ielem,:) = [inode   inode+1   inode+(npx)];
                T(ielem+1,:) = [inode+1+npx   inode+npx inode+1];
                
            end   
        end
        end
    end
       
        
end
%Noeuds sur le bord stockées dans "b" comme suit (bas,gauche,droite,haut)
%P1
if n==0
X1=zeros((2*npx-1)*(npy-1)+(nx*npy),2);T1=zeros(nx*ny,3);
b=[1:npx,npx+1:npx:npx*npy,2*npx:npx:npx*npy,npx*npy-npx+2:npx*npy-1];
trimesh(T,X(:,1),X(:,2));
hold on
x=X(:,1);y=X(:,2);
plot(x,y,'o');
end
if n==1
    if k==1
%NP1
b=[1:nx,nx+1:2*npx-1+nx:nx*npy+ny*(2*npx-1),2*npx-1+nx:2*npx-1+nx:nx*npy+ny*(2*npx-1),(2*npx-1+nx)*ny+1:(2*npx-1)*(npy-1)+(nx*npy)];
trimesh(T,X(:,1),X(:,2));
hold on
x=X1(:,1);y=X1(:,2);
plot(x,y,'o');    
    elseif k==2
%NP2    
b=[1:2*nx,sort([2*nx+1:2*(2*npx-1)+2*nx:2*(2*npx-1)*ny+(2*nx*npy),(2*npx-1)+2*nx+1:2*(2*npx-1)+2*nx:2*(2*npx-1)*ny+(2*nx*npy)]),sort([2*nx+(2*npx-1):2*(2*npx-1)+2*nx:2*(2*npx-1)*ny+(2*nx*npy),2*(2*npx-1)+2*nx:2*(2*npx-1)+2*nx:2*(2*npx-1)*ny+(2*nx*npy)]),sort(length(X):-1:length(X)-2*nx+1)];
trimesh(T,X(:,1),X(:,2));
hold on
x=X1(:,1);y=X1(:,2);
plot(x,y,'o');
    end
%
%
%}
end
end


