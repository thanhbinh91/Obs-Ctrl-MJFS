function [gama,Lm,Fm,Sm,Zm] = LMIs_Ex1(beta,alpha)
[nx, nu, nw, ny, nz, r, A, B, E, C, D, G, H, J] = SysParas_Ex1;

nlmis = 1;
setlmis([]);

ga2 = lmivar(1,[1 1]);
lmiterm([-nlmis,   1,   1,   ga2],1,1);
nlmis = nlmis + 1;

P1 = zeros(r,r);
P2 = zeros(r,r);
P3 = zeros(r,r);


S  = zeros(1,r);
X  = zeros(1,r);
Y  = zeros(1,r);
Z  = zeros(1,r);
F_ = zeros(1,r);
L_ = zeros(1,r);
W  = zeros(1,r);
U  = zeros(r,r);

Q1 = lmivar(1, [nx,  1]);
Q2 = lmivar(1, [nx,  1]);


for i = 1:r
    X(i)  = lmivar(2, [nz,  nz]);
    Y(i)  = lmivar(2, [nx,  nx]);
    Z(i)  = lmivar(2, [nx,  nx]);
    F_(i) = lmivar(2, [nu,  nx]);
    L_(i) = lmivar(2, [nx,  ny]);
    S(i)  = lmivar(2, [nu,  nu]);
    W(i)  = lmivar(2, [2*nx+nz+nu, 4*nx+nw+nz+nu]);
    for j = 1:r
        P1(i,j) = lmivar(1, [nx,  1]);
        P2(i,j) = lmivar(1, [nx,  1]);
        P3(i,j) = lmivar(2, [nx,  nx]);
        U(i,j)  = lmivar(1, [2*nx+nz+nu ,1]);
    end
end



for i = 1:r
    lmiterm([nlmis  1  1  Q1],      -1 ,1);
    lmiterm([nlmis  1  1  P1(i,i)],  1 ,1);    
    lmiterm([nlmis  2  2  Q2],      -1 ,1);
    lmiterm([nlmis  2  2  P2(i,i)],  1 ,1);
    lmiterm([nlmis  2  1  P3(i,i)],  1 ,1); 
    nlmis = nlmis + 1;
    
    Up(1,i,i);
    nlmis = nlmis + 1;
for j = 1:r
    if j ~= i
    lmiterm([nlmis  1  1  Q1],      -r/(r-1) ,1);
    lmiterm([nlmis  1  1  P1(i,i)],  1/(r-1) ,1);
    lmiterm([nlmis  1  1  P1(i,j)],  1/2     ,1); 
    lmiterm([nlmis  1  1  P1(j,i)],  1/2     ,1); 
    lmiterm([nlmis  2  2  Q2],      -r/(r-1) ,1);
    lmiterm([nlmis  2  2  P2(i,i)],  1/(r-1) ,1);
    lmiterm([nlmis  2  2  P2(i,j)],  1/2     ,1); 
    lmiterm([nlmis  2  2  P2(j,i)],  1/2     ,1); 
    lmiterm([nlmis  2  1  P3(i,i)],  1/(r-1) ,1);
    lmiterm([nlmis  2  1  P3(i,j)],  1/2     ,1); 
    lmiterm([nlmis  2  1  P3(j,i)],  1/2     ,1); 
    nlmis = nlmis + 1;
    
    Up(1/(r-1),i,i);
    Up(1/2,i,j);
    Up(1/2,j,i);
    nlmis = nlmis + 1;
    end
end
end


%%%%% LMI solver %%%%%
lmis = getlmis;
options = [1e-3,300,0,0,0];
c = zeros(1,decnbr(lmis)); c(1) = 1;
[copt,xopt] = mincx(lmis,c,options);
% [t,xopt] = feasp(lmis,options);
gama = sqrt(dec2mat(lmis, xopt, ga2));

Lm = zeros(nx,ny,r);
Fm = zeros(nu,nx,r);
Sm = zeros(nu,nu,r);
Zm = zeros(nx,nx,r);
Um = zeros(6,6,r,r);
P1m = zeros(nx,nx,r,r);
P2m = zeros(nx,nx,r,r);
Xm = zeros(nz,nz,r);
Ym = zeros(nx,nx,r);
Wm = zeros(2*nx+nz+nu, 4*nx+nw+nz+nu,r);

for i = 1:r
    Lm(:,:,i) = dec2mat(lmis, xopt, L_(i));
    Fm(:,:,i) = dec2mat(lmis, xopt, F_(i));
    Sm(:,:,i) = dec2mat(lmis, xopt, S(i));
    Zm(:,:,i) = dec2mat(lmis, xopt, Z(i));
    Xm(:,:,i) = dec2mat(lmis, xopt, X(i));
    Ym(:,:,i) = dec2mat(lmis, xopt, Y(i));
    Wm(:,:,i) = dec2mat(lmis, xopt, W(i));
    for j = 1:r
        Um(:,:,i,j)  = dec2mat(lmis, xopt, U(i,j));
        P1m(:,:,i,j) = dec2mat(lmis, xopt, P1(i,j));
        P2m(:,:,i,j) = dec2mat(lmis, xopt, P2(i,j));
    end
end


function Up(co,i,j)
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
    lmiterm([nlmis    1   1    U(l,i)], -alpha^2*Om1', Om1*co);
    end
    %% _Up(1,1)
    % Up1
    lmiterm([nlmis    1   1    ga2],     -e3',    e3*co);
    lmiterm([nlmis    1   1    0],        e4'*e4*co);
    lmiterm([nlmis    1   1    Q1],       e5',    e5*co);
    lmiterm([nlmis    1   1    Q2],       e6',    e6*co);
    

    % Up2
    lmiterm([nlmis    1   1    X(i)],    -e4',    e4*co, 's');
    lmiterm([nlmis    1   1    Y(i)],    -e5',    e5*co, 's');
    lmiterm([nlmis    1   1    Z(i)],    -e6',    e6*co, 's');
    lmiterm([nlmis    1   1    F_(i)],    e7',    e1*co, 's');
    lmiterm([nlmis    1   1    F_(i)],   -e7',    e1*co, 's');
    

    % Up4
    lmiterm([nlmis    1   1    P1(i,j)], -e1',    e1*co);
    lmiterm([nlmis    1   1    P2(i,j)], -e2',    e2*co);
    lmiterm([nlmis    1   1    P3(i,j)], -e2',    e1*co,'s');
    lmiterm([nlmis    1   1    Z(j)],    -e6',    A(:,:,i)*e1*co,'s');
    lmiterm([nlmis    1   1    L_(j)],   -e6',    C(:,:,i)*e1*co,'s');
    lmiterm([nlmis    1   1    Z(j)],     e6',    A(:,:,i)*e2*co,'s');
    lmiterm([nlmis    1   1    L_(j)],   -e6',    C(:,:,i)*e2*co,'s');
    lmiterm([nlmis    1   1   -Z(j)],  -beta*e7'*B(:,:,i)',  e6*co,'s');
    lmiterm([nlmis    1   1    S(j)],  -beta*e7',            e7*co,'s');


    % Up3
    lmiterm([nlmis    1   1    X(j)],     e4',                  G(:,:,i)*e1*co,'s');
    lmiterm([nlmis    1   1    F_(j)],    e4'*H(:,:,i),         e1*co,'s');
    lmiterm([nlmis    1   1    F_(j)],   -e4'*H(:,:,i),         e2*co,'s');
    lmiterm([nlmis    1   1    X(j)],     e4',                  J(:,:,i)*e3*co,'s'); 

    lmiterm([nlmis    1   1    Y(j)],     e5',                  A(:,:,i)*e1*co,'s');
    lmiterm([nlmis    1   1    F_(j)],    e5'*B(:,:,i),         e1*co,'s');
    lmiterm([nlmis    1   1    F_(j)],   -e5'*B(:,:,i),         e2*co,'s');
    lmiterm([nlmis    1   1    Y(j)],     e5',                  E(:,:,i)*e3*co,'s');

    lmiterm([nlmis    1   1    Z(j)],     e6',                  A(:,:,i)*e1*co,'s');
    lmiterm([nlmis    1   1    L_(j)],    e6',                  C(:,:,i)*e1*co,'s');
    lmiterm([nlmis    1   1    Z(j)],     e6',                  E(:,:,i)*e3*co,'s');
    lmiterm([nlmis    1   1    L_(j)],   -e6',                  D(:,:,i)*e3*co,'s');

    lmiterm([nlmis    1   1   -X(j)],     beta*e7'*H(:,:,i)',   e4*co,'s');
    lmiterm([nlmis    1   1   -S(j)],    -beta*e7',             H(:,:,i)'*e4*co,'s');
    lmiterm([nlmis    1   1   -Y(j)],     beta*e7'*B(:,:,i)',   e5*co,'s');
    lmiterm([nlmis    1   1   -S(j)],    -beta*e7',             B(:,:,i)'*e5*co,'s');
    lmiterm([nlmis    1   1   -Z(j)],     beta*e7'*B(:,:,i)',   e6*co,'s');

    %% _Up(2,1) & _Up(2,2)
    c1 = [eye(nz)           zeros(nz,2*nx+nu)];
    c2 = [zeros(nx, nz)     eye(nx) zeros(nx,nx+nu)];
    c3 = [zeros(nx, nz+nx)  eye(nx) zeros(nx,nu)];
    c4 = [zeros(nu,nz+2*nx) eye(nu)];
    

    for l = 1:r
    lmiterm([nlmis   1+l  1    X(i)],     c1',                  G(:,:,l)*e1*co);
    lmiterm([nlmis   1+l  1    F_(i)],    c1'*H(:,:,l),         e1*co);
    lmiterm([nlmis   1+l  1    F_(i)],   -c1'*H(:,:,l),         e2*co);
    lmiterm([nlmis   1+l  1    X(i)],     c1',                  J(:,:,l)*e3*co); 

    lmiterm([nlmis   1+l  1    Y(i)],     c2',                  A(:,:,l)*e1*co);
    lmiterm([nlmis   1+l  1    F_(i)],    c2'*B(:,:,l),         e1*co);
    lmiterm([nlmis   1+l  1    F_(i)],   -c2'*B(:,:,l),         e2*co);
    lmiterm([nlmis   1+l  1    Y(i)],     c2',                  E(:,:,l)*e3*co);

    lmiterm([nlmis   1+l  1    Z(i)],     c3',                  A(:,:,l)*e1*co);
    lmiterm([nlmis   1+l  1    L_(i)],    c3',                  C(:,:,l)*e1*co);
    lmiterm([nlmis   1+l  1    Z(i)],     c3',                  E(:,:,l)*e3*co);
    lmiterm([nlmis   1+l  1    L_(i)],   -c3',                  D(:,:,l)*e3*co);

    lmiterm([nlmis   1+l  1   -X(i)],     beta*c4'*H(:,:,l)',   e4*co);
    lmiterm([nlmis   1+l  1   -S(i)],    -beta*c4',             H(:,:,l)'*e4*co);
    lmiterm([nlmis   1+l  1   -Y(i)],     beta*c4'*B(:,:,l)',   e5*co);
    lmiterm([nlmis   1+l  1   -S(i)],    -beta*c4',             B(:,:,l)'*e5*co);
    lmiterm([nlmis   1+l  1   -Z(i)],     beta*c4'*B(:,:,l)',   e6*co);

    lmiterm([nlmis   1+l  1    W(i)],     1,                    co);

    lmiterm([nlmis   1+l 1+l   U(l,i)],   1,                    co);
    end
 
end

end