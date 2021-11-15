function [gama,Lm,Fm,Sm,Zm] = LMIs_Ex2(beta,alpha)
[nx, nu, nw, ny, nz, s, r, A, B, E, C, D, G, H, J, Hk, Hb, Up_Pi, Lo_Pi] = SysParas_Ex2;

nlmis = 1;
setlmis([]);

ga2 = lmivar(1,[1 1]);
lmiterm([-nlmis,   1,   1,   ga2], 1, 1);
nlmis = nlmis + 1;

P1 = zeros(s,r,r);
P2 = zeros(s,r,r);
P3 = zeros(s,r,r);


S  = zeros(s,r);
X  = zeros(s,r);
Y  = zeros(s,r);
Z  = zeros(s,r);
F_ = zeros(s,r);
L_ = zeros(s,r);
W1 = zeros(s,r);
W2 = zeros(s,r);
U  = zeros(s,r,r);

Q1 = zeros(1,s);
Q2 = zeros(1,s);
T = zeros(s,s);




for g = 1:s
    Q1(g) = lmivar(1, [nx,  1]);
    Q2(g) = lmivar(1, [nx,  1]);

    for hh = 1:s
        if hh~= g
            T(g,hh) = lmivar(2, [2*nx,  2*nx]);
        end
    end

    for i = 1:r
        X(g,i)  = lmivar(2, [nz,  nz]);
        Y(g,i)  = lmivar(2, [nx,  nx]);
        Z(g,i)  = lmivar(2, [nx,  nx]);
        F_(g,i) = lmivar(2, [nu,  nx]);
        L_(g,i) = lmivar(2, [nx,  ny]);
        S(g,i)  = lmivar(2, [nu,  nu]);
        W1(g,i) = lmivar(2, [2*nx+nz, 2*nx+nw]);
        W2(g,i) = lmivar(2, [nu, 2*nx+nz]);
        for j = 1:r
            P1(g,i,j) = lmivar(1, [nx,  1]);
            P2(g,i,j) = lmivar(1, [nx,  1]);
            P3(g,i,j) = lmivar(2, [nx,  nx]);
            U(g,i,j)  = lmivar(1, [2*nx+nz+nu ,1]);
        end
    end
end

%%
for g = 1:s
    for i = 1:r
        Phi(1,g,i,i);
        nlmis = nlmis + 1;

        Up(1,g,i,i);
        nlmis = nlmis + 1;
    for j = 1:r
        if j ~= i

        Phi(1/(r-1),g,i,i);
        Phi(1/2,g,i,j);
        Phi(1/2,g,j,i);
        nlmis = nlmis + 1;
        
        
        Up(1/(r-1),g,i,i);
        Up(1/2,g,i,j);
        Up(1/2,g,j,i);
        nlmis = nlmis + 1;
        end
    end
    end
end
%%


    %%%%% LMI solver %%%%%
    lmis = getlmis;
    options = [1e-2,1000,0,0,0];
    c = zeros(1,decnbr(lmis)); c(1) = 1;
    [copt,xopt] = mincx(lmis,c,options);
    % [t,xopt] = feasp(lmis,options);
    gama = sqrt(dec2mat(lmis, xopt, ga2));
    
    Lm = zeros(nx,ny,s,r);
    Fm = zeros(nu,nx,s,r);
    Sm = zeros(nu,nu,s,r);
    Zm = zeros(nx,nx,s,r);
    for g = 1:s
    for i = 1:r
        Lm(:,:,g,i) = dec2mat(lmis, xopt, L_(g,i));
        Fm(:,:,g,i) = dec2mat(lmis, xopt, F_(g,i));
        Sm(:,:,g,i) = dec2mat(lmis, xopt, S(g,i));
        Zm(:,:,g,i) = dec2mat(lmis, xopt, Z(g,i));
    end
    end




    function Up(co,g,i,j)
    %     B_  =  B;
        Om1 = [zeros(2*nx+nz+nu, 2*nx+nw)  eye(2*nx+nz+nu)];
        e1  = [eye(nx)              zeros(nx,3*nx+nz+nu+nw)];
        e2  = [zeros(nx)            eye(nx) zeros(nx, 2*nx+nz+nu+nw)];
        e3  = [zeros(nw,2*nx)       eye(nw) zeros(nw, 2*nx+nz+nu)];
        e4  = [zeros(nz,2*nx+nw)    eye(nz) zeros(nz, 2*nx+nu)];
        e5  = [zeros(nx,2*nx+nw+nz) eye(nx) zeros(nx, nx+nu)];
        e6  = [zeros(nx,3*nx+nw+nz) eye(nx) zeros(nx, nu)];
        e7  = [zeros(nu,4*nx+nw+nz) eye(nu)]; 
    
        for l = 1:r
        lmiterm([nlmis    1   1    U(g,l,i)], -alpha^2*Om1', Om1*co);
        end
        %% _Up(1,1)
        % Up1
        lmiterm([nlmis    1   1    ga2],       -e3',    e3*co);
        lmiterm([nlmis    1   1    0],          e4'*e4*co);
        lmiterm([nlmis    1   1    Q1(g)],      e5',    e5*co);
        lmiterm([nlmis    1   1    Q2(g)],      e6',    e6*co);
        
    
        % Up2
        lmiterm([nlmis    1   1    X(g,i)],    -e4',    e4*co, 's');
        lmiterm([nlmis    1   1    Y(g,i)],    -e5',    e5*co, 's');
        lmiterm([nlmis    1   1    Z(g,i)],    -e6',    e6*co, 's');
        lmiterm([nlmis    1   1    F_(g,i)],    e7',    e1*co, 's');
        lmiterm([nlmis    1   1    F_(g,i)],   -e7',    e1*co, 's');
        
    
        % Up4
        lmiterm([nlmis    1   1    P1(g,i,j)], -e1',    e1*co);
        lmiterm([nlmis    1   1    P2(g,i,j)], -e2',    e2*co);
        lmiterm([nlmis    1   1    P3(g,i,j)], -e2',    e1*co, 's');
        lmiterm([nlmis    1   1    Z(g,j)],    -e6',    A(:,:,g,i)*e1*co,'s');
        lmiterm([nlmis    1   1    L_(g,j)],   -e6',    C(:,:,g,i)*e1*co,'s');
        lmiterm([nlmis    1   1    Z(g,j)],     e6',    A(:,:,g,i)*e2*co,'s');
        lmiterm([nlmis    1   1    L_(g,j)],   -e6',    C(:,:,g,i)*e2*co,'s');
        lmiterm([nlmis    1   1   -Z(g,j)],    -beta*e7'*B(:,:,g,i)',  e6*co,'s');
        lmiterm([nlmis    1   1    S(g,j)],    -beta*e7',              e7*co,'s');
    
    
        % Up3
        lmiterm([nlmis    1   1    X(g,j)],     e4',                  G(:,:,g,i)*e1*co,'s');
        lmiterm([nlmis    1   1    F_(g,j)],    e4'*H(:,:,g,i),       e1*co,'s');
        lmiterm([nlmis    1   1    F_(g,j)],   -e4'*H(:,:,g,i),       e2*co,'s');
        lmiterm([nlmis    1   1    X(g,j)],     e4',                  J(:,:,g,i)*e3*co,'s'); 
    
        lmiterm([nlmis    1   1    Y(g,j)],     e5',                  A(:,:,g,i)*e1*co,'s');
        lmiterm([nlmis    1   1    F_(g,j)],    e5'*B(:,:,g,i),       e1*co,'s');
        lmiterm([nlmis    1   1    F_(g,j)],   -e5'*B(:,:,g,i),       e2*co,'s');
        lmiterm([nlmis    1   1    Y(g,j)],     e5',                  E(:,:,g,i)*e3*co,'s');
    
        lmiterm([nlmis    1   1    Z(g,j)],     e6',                  A(:,:,g,i)*e1*co,'s');
        lmiterm([nlmis    1   1    L_(g,j)],    e6',                  C(:,:,g,i)*e1*co,'s');
        lmiterm([nlmis    1   1    Z(g,j)],     e6',                  E(:,:,g,i)*e3*co,'s');
        lmiterm([nlmis    1   1    L_(g,j)],   -e6',                  D(:,:,g,i)*e3*co,'s');
    
        lmiterm([nlmis    1   1   -X(g,j)],     beta*e7'*H(:,:,g,i)', e4*co,'s');
        lmiterm([nlmis    1   1   -S(g,j)],    -beta*e7',             H(:,:,g,i)'*e4*co,'s');
        lmiterm([nlmis    1   1   -Y(g,j)],     beta*e7'*B(:,:,g,i)', e5*co,'s');
        lmiterm([nlmis    1   1   -S(g,j)],    -beta*e7',             B(:,:,g,i)'*e5*co,'s');
        lmiterm([nlmis    1   1   -Z(g,j)],     beta*e7'*B(:,:,g,i)', e6*co,'s');
    
        %% _Up(2,1) & _Up(2,2)
        c1 = [eye(nz)           zeros(nz,2*nx+nu)];
        c2 = [zeros(nx, nz)     eye(nx) zeros(nx,nx+nu)];
        c3 = [zeros(nx, nz+nx)  eye(nx) zeros(nx,nu)];
        c4 = [zeros(nu,nz+2*nx) eye(nu)];
        
    
        for l = 1:r
        lmiterm([nlmis   1+l  1    X(g,i)],     c1',                  G(:,:,g,l)*e1*co);
        lmiterm([nlmis   1+l  1    F_(g,i)],    c1'*H(:,:,g,l),       e1*co);
        lmiterm([nlmis   1+l  1    F_(g,i)],   -c1'*H(:,:,g,l),       e2*co);
        lmiterm([nlmis   1+l  1    X(g,i)],     c1',                  J(:,:,g,l)*e3*co); 
    
        lmiterm([nlmis   1+l  1    Y(g,i)],     c2',                  A(:,:,g,l)*e1*co);
        lmiterm([nlmis   1+l  1    F_(g,i)],    c2'*B(:,:,g,l),       e1*co);
        lmiterm([nlmis   1+l  1    F_(g,i)],   -c2'*B(:,:,g,l),       e2*co);
        lmiterm([nlmis   1+l  1    Y(g,i)],     c2',                  E(:,:,g,l)*e3*co);
    
        lmiterm([nlmis   1+l  1    Z(g,i)],     c3',                  A(:,:,g,l)*e1*co);
        lmiterm([nlmis   1+l  1    L_(g,i)],    c3',                  C(:,:,g,l)*e1*co);
        lmiterm([nlmis   1+l  1    Z(g,i)],     c3',                  E(:,:,g,l)*e3*co);
        lmiterm([nlmis   1+l  1    L_(g,i)],   -c3',                  D(:,:,g,l)*e3*co);
    
        lmiterm([nlmis   1+l  1   -X(g,i)],     beta*c4'*H(:,:,g,l)', e4*co);
        lmiterm([nlmis   1+l  1   -S(g,i)],    -beta*c4',             H(:,:,g,l)'*e4*co);
        lmiterm([nlmis   1+l  1   -Y(g,i)],     beta*c4'*B(:,:,g,l)', e5*co);
        lmiterm([nlmis   1+l  1   -S(g,i)],    -beta*c4',             B(:,:,g,l)'*e5*co);
        lmiterm([nlmis   1+l  1   -Z(g,i)],     beta*c4'*B(:,:,g,l)', e6*co);
        

        q1 = [eye(2*nx+nz)           zeros(2*nx+nz,nu)];
        q2 = [zeros(nu,2*nx+nz),     eye(nu)];
        f1 = [eye(2*nx+nw)           zeros(2*nx+nw,2*nx+nz+nu)];
        f2 = [zeros(2*nx+nz,2*nx+nw+nu) eye(2*nx+nz)];


        lmiterm([nlmis   1+l  1    W1(g,i)],     q1',                 f1*co);
        lmiterm([nlmis   1+l  1    W2(g,i)],     q2',                 f2*co);
    
        lmiterm([nlmis   1+l 1+l   U(g,l,i)],   1,                    co);
        end
     
    end



    function Phi(co,g,i,j)
        e1 = [eye(nx)   zeros(nx)];
        e2 = [zeros(nx) eye(nx)];

        lmiterm([nlmis  1  1  Q1(g)],       -e1',            e1*co);
        lmiterm([nlmis  1  1  P1(g,i,j)],    e1',            e1*co);
        lmiterm([nlmis  1  1  P3(g,i,j)],    e2',            e1*co, 's');
        
        lmiterm([nlmis  1  1  Q2(g)],       -e2',            e2*co);
        lmiterm([nlmis  1  1  P2(g,i,j)],    e2',            e2*co);

        for h = setdiff(Hk{g},g,'legacy')
            lmiterm([nlmis  1  1  P1(h,i,j)],    e1'*Up_Pi(g,h), e1*co);
            lmiterm([nlmis  1  1  P2(h,i,j)],    e2'*Up_Pi(g,h), e2*co);
            lmiterm([nlmis  1  1  P3(h,i,j)],    e2'*Up_Pi(g,h), e1*co,'s');
            lmiterm([nlmis  1  1  P1(g,i,j)],   -e1'*Up_Pi(g,h), e1*co);
            lmiterm([nlmis  1  1  P2(g,i,j)],   -e2'*Up_Pi(g,h), e2*co);
            lmiterm([nlmis  1  1  P3(g,i,j)],   -e2'*Up_Pi(g,h), e1*co,'s');
        end

        for h = setdiff(Hb{g},g,'legacy')
            lmiterm([nlmis  1  1  T(g,h)],   Up_Pi(g,h)*Lo_Pi(g,h), co,'s');
        end
        
        row = 0;
        for h = setdiff(Hb{g},g,'legacy')
            row = row+1;
            lmiterm([nlmis  1+row  1  P1(h,i,j)],   e1',            e1*co);
            lmiterm([nlmis  1+row  1  P2(h,i,j)],   e2',            e2*co);
            lmiterm([nlmis  1+row  1  P3(h,i,j)],   e2',            e1*co);
            lmiterm([nlmis  1+row  1 -P3(h,i,j)],   e1',            e2*co);

            lmiterm([nlmis  1+row  1  P1(g,i,j)],  -e1',            e1*co);
            lmiterm([nlmis  1+row  1  P2(g,i,j)],  -e2',            e2*co);
            lmiterm([nlmis  1+row  1  P3(g,i,j)],  -e2',            e1*co);
            lmiterm([nlmis  1+row  1 -P3(g,i,j)],  -e1',            e2*co);


            lmiterm([nlmis  1+row  1      T(g,h)], -Up_Pi(g,h)-Lo_Pi(g,h),  co);
            lmiterm([nlmis  1+row  1+row  T(g,h)], 1,                      co,'s');
        end
    end

end
