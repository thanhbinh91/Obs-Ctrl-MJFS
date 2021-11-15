function [nx, nu, nw, ny, nz, s, r, A, B, E, C, D, G, H, J, Hk, Hb, Up_Pi, Lo_Pi] = SysParas_Ex2
    %% Dimension 
    nx = 2;
    nu = 1; 
    nw = 1; 
    ny = 1; 
    nz = 1; 
    r  = 2;
    s  = 3;


    %% MJFS
    A  = zeros(nx,nx,s,r);
    B  = zeros(nx,nu,s,r);
    E  = zeros(nx,nw,s,r);
    C  = zeros(ny,nx,s,r);
    D  = zeros(ny,nw,s,r);
    G  = zeros(nz,nx,s,r);
    H  = zeros(nz,nu,s,r);
    J  = zeros(nz,nw,s,r);


    Hk = cell(1,3);
    Hb = cell(1,3);
    Hk{1} = [1 2 3]; Hk{2} = 3;     Hk{3} = 3;
    Hb{1} = [];      Hb{2} = [1 2]; Hb{3} = [1 2];

    Pi    = [.8  .1  .1;
             .3  .6  .1;
             .5  .2  .3];
    Up_Pi = [.8  .1  .1;
             .5  .8  .1;
             .65 .35 .3];
    Lo_Pi = [.8  .1  .1;
             .1  .4  .1;
             .35 .05 .3];
    
    a = [-2 -1.5 -1];
    b = [-1  0    1];

    %% Markov chain desription   
    for g = 1:s
        A(:,:,g,1) =  [1+a(g) -0.5;
                       1       0];
              
        A(:,:,g,2) =  [-1     -0.5;
                        1      0];
              
        B(:,:,g,1) = [ 1;  1-b(g)];
        B(:,:,g,2) = [-2;  1];
    
        E(:,:,g,1) = [0.2; 0.3];
        E(:,:,g,2) = [0.5;-0.1];
    
        C(:,:,g,1) = [0.1 -0.4];           
        C(:,:,g,2) = [-0.2 -0.6];
    
        D(:,:,g,1) = -0.1;
        D(:,:,g,2) = -0.2;
    
        G(:,:,g,1) = [1   0.5];
        G(:,:,g,2) = [0.5 1];
    
        H(:,:,g,1) = 1;
        H(:,:,g,2) = 0.5;

        J(:,:,g,1) = 0.4;
        J(:,:,g,2) = 0.2;
    end

end