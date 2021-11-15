function [nx, nu, nw, ny, nz, r, A, B, E, C, D, G, H, J] = SysParas_Ex1
%% Dimension 
nx = 2;
nu = 1; 
nw = 1; 
ny = 1; 
nz = 1; 
r  = 2;
% s  = 1; %No Markov

%% MJFS
A  = zeros(nx,nx,r);
B  = zeros(nx,nu,r);
E  = zeros(nx,nw,r);
C  = zeros(ny,nx,r);
D  = zeros(ny,nw,r);
G  = zeros(nz,nx,r);
H  = zeros(nz,nu,r);
J  = zeros(nz,nw,r);




%% Markov chain desription
%% 
    a = -1.6;
    b = 0;

    A(:,:,1) = [1+a -0.5;
                1    0];
              
    A(:,:,2) = [-1  -0.5;
                 1   0];
              
    B(:,:,1) = [ 1;  1-b];
    B(:,:,2) = [-2;  1];
    
    E(:,:,1) = [0.2; 0.3];
    E(:,:,2) = [0.5;-0.1];
    
    C(:,:,1) = [0.1 -0.4];         
    C(:,:,2) = [-0.2 -0.6];
    
    D(:,:,1) = -0.1;
    D(:,:,2) = -0.2;
    
    G(:,:,1)  = [ 1   0.5];
    G(:,:,2)  = [ 0.5 1];
    
    H(:,:,1)  = 1;
    H(:,:,2)  = 0.5;
     
    J(:,:,1) =  0.4;
    J(:,:,2) =  0.2;
   

end