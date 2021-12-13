% 2 dimensional particle state space size (raduis) and density.
% 
N = 30 % size classes enumerated [0, 1, ..., N-1]
M = 10 % density classes enumerated [0, 1, ..., M-1]
L = N*M % number of bins: b index [0, 1, ..., L-1]
K = (L+1)*L/2 % number of combos: k index [0, 1, ..., K-1]
k = [0:K-1]';
b = [0:L-1]';
z = (2*L + 1 - sqrt((2*L + 1).^2 - 8*k))/2;
bi = floor(z);
bj = k - bi*L + bi.*(bi - 1)/2 + bi;
xi = floor(bi/M); zi = bi - xi*M;
xj = floor(bj/M); zj = bj - xj*M;
x = [0:N-1]; z = [0:M-1];

ro = 1;             % smallest size [µ]
a = 2.7;            % self-similar parameter; rioj^a = ri^a + rj^a
%delta = 2^(1/a);    % size class incremant factor; ri = ro delta^i
delta = 1.10;        % natural fro self similarity; delta = 2^(1/a) 
dt = 1;            % time step [day]
TotalP = 1          % total productivity [g /m^2 / day] (carbon and ballast)
q = delta^(a-3);
nu = 1E-6           % viscosity of water [m2/s]
rho_water = 1E-6;   %density of water [µg/µm^3];
d_rho = 0.4*rho_water/(M + 1)
pip = 4*pi/3;

xz = @(b) [floor(b/(M)), b - floor(b/(M))*(M)];
logd = @(x) log(x)./log(delta);
zeta = @(xi,xj) (1 + delta.^(a*(xi - xj)));
p = @(x,z) x*(M) + z;


rad = @(x) ro*delta.^x;
den = @(x,z) d_rho*z.*q.^x;
mass = @(x,z) pip*(den(x,z) + rho_water).*(rad(x).^3);
wsink = @(x,z) 2/9*9.8/nu*den(x,z)/rho_water .* ((rad(x)*1E-6).^2) * 24 * 3600 % sinking speed [m / day]
%Re = w.*(r*1E-6)/nu; % Reynolds number
beta = 1E-5; % constant size independent encounter rate
H = 50; % depth of the surface mixed layer [m]

xioj = xi + logd(zeta(xj,xi))/a;
zioj = zi./zeta(xj,xi) + zj./zeta(xi,xj);

x300 = floor(xioj);
z300 = floor(zioj);
b300 = p(x300,z300);
b310 = p(x300+1,z300);
b301 = p(x300,z300+1);
b311 = p(x300+1,z300+1);
dx1 = xioj - x300;
dx0 = 1 - dx1;
dz1 = zioj - z300;
dz0 = 1 - dz1;
f00 = dx0.*dz0;
f10 = dx1.*dz0;
f01 = dx0.*dz1;
f11 = dx1.*dz1;

xzs = xz(b);
m = mass(xzs(:,1), xzs(:,2));  m = reshape(m,M,N);
w = wsink(xzs(:,1), xzs(:,2)); w = reshape(w,M,N);
r = rad(xzs(:,1));
beta = 1E-7;
betaw = beta*ones(1,K);
betat = beta*ones(1,K);
%calculate coagulation rates
for k = 1:K
    % differential settling
    b1 = bi(k)+1;   b2 = bj(k)+1;
    m1 = m(b1);     m2 = m(b2);
    w1 = w(b1);     w2 = w(b2);
    x1 = xi(k);     x2 = xj(k);
    r1 = rad(x1)*1E-6;   r2 = rad(x2)*1E-6;
    betaw(k) = 0.5*pi*min(r1,r2)^2*abs(w1-w2); % m^3/day
    epsilon = 1E-4 % turbulent dissipation rate [m^2 / s^3]
    betat(k) = 24*3600*9*sqrt(epsilon/nu)*(r1 + r2)^3; % m^3/day
end

Numb = zeros(M,N);
dNum = zeros(M,N);
Mass = zeros(M,N);
dMass = zeros(M,N);
Prod = zeros(M,N);
Prod(1:M-1,1) = 0.1/H;
Mass = Mass + Prod;
Numb = Mass./m
figure(1)
clf;
subplot(3,1,1)
matplot1 = imagesc(x,z,Mass);title('Mass')
timtext = text(.8,1.1,['time ' int2str(0)],'Units','normalized')
axis xy;
subplot(3,1,2)
matplot2 = imagesc(x,z,Numb);title('Number')
axis xy;
subplot(3,1,3)
matplot3 = plot(x+.5,zeros(size(x)),'o-');title('Flux');
fluxtext = text(.8,1.1,['flux ' num2str(0,'%5.3f')],'Units','normalized')


for t = 1:2000
     Numb = Mass./m;
     dMass(:,:) = 0;
     for k = 1:K % for each possible combination of boxes
        b1 = bi(k)+1; % decompose the two indices of the two boxes involved
        b2 = bj(k)+1;
        dN = beta*Numb(b1).*Numb(b2); % the number of particles formed 
        if dN > 0; % only bother if dN >0
            m1 = m(b1);
            m2 = m(b2);
            bk00 = b300(k) + 1;
            if bk00 > L; bk00 = L; end;
            bk01 = bk00 + 1;
            if bk01 > L; bk01 = bk00; end;
            bk10 = bk00 + M;
            if bk10 > L; bk10 = bk00; end;
            bk11 = bk01 + M;
            if bk11 > L; bk11 = bk01; end;
            %dNum
            dMass(b1) = dMass(b1) - dN*m1;
            dMass(b2) = dMass(b2) - dN*m2;
            dMass(bk00) = dMass(bk00) + f00(k)*dN*(m1 + m2);
            dMass(bk01) = dMass(bk01) + f01(k)*dN*(m1 + m2);
            dMass(bk10) = dMass(bk10) + f10(k)*dN*(m1 + m2);
            dMass(bk11) = dMass(bk11) + f11(k)*dN*(m1 + m2);
        end
     end
     sum(sum(dMass));
     Mass = Mass + dt*(dMass + Prod);
     Mass(Mass<0)=0;
     dMass_w = (w*dt/H).*Mass;
     q = sum(dMass_w);
     Q = sum(q);
     Mass = Mass - dMass_w;
     set(matplot1,'Cdata',Mass)
     set(matplot2,'Cdata',Numb)
     set(matplot3,'ydata',q)

     if mod(t,10)==0
        set(timtext,'String',['time ' int2str(t)]);
        set(fluxtext,'String',['flux ' num2str(Q,'%6.4f')]);
     end 
     drawnow
end
