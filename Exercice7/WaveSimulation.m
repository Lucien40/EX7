%% Wave simulation and analysis
%% Parameter initialisation
repertoire = './';
executable = 'Exercice7';
input = 'configurationT.in';
param = 'schema';
paramval = ['A' 'B' 'C'];
paramnum = size(paramval,2);

%% Tsunami A, B, C simulations 
for i = 1:paramnum
    cmd = sprintf('wsl %s%s %s %s=%s %s=%s', repertoire, executable, input,...
                                             param, paramval(i), 'output', [param '=' paramval(i)]);
    disp(cmd); system(cmd);
end

%% Load data
fichier = 'schema=A';
data = load([fichier,'_u.out']);
x = data(:,1);
u = data(:,2);
data = load([fichier,'_E.out']);
t = data(:,1);
E = data(:,2);
data = load([fichier,'_f.out']);
f = data(:,2:end);
fichier = 'schema=B';
data = load([fichier,'_u.out']);
x2 = data(:,1);
u2 = data(:,2);
data = load([fichier,'_E.out']);
t2 = data(:,1);
E2 = data(:,2);
data = load([fichier,'_f.out']);
f2 = data(:,2:end);
fichier = 'schema=C';
data = load([fichier,'_u.out']);
x3 = data(:,1);
u3 = data(:,2);
data = load([fichier,'_E.out']);
t3 = data(:,1);
E3 = data(:,2);
data = load([fichier,'_f.out']);
f3 = data(:,2:end);

h_ocean = 8e3;
h_recif = 20;
xa = 200e3;
xb = 370e3;
xc = 430e3;
xd = 600e3;
g = 9.81;
h = zeros(size(x,1),1);
h(x <= xa | x >= xd) = h_ocean;
h(xa < x & x < xb) = h_ocean + (h_recif - h_ocean) * ...
                     sin(0.5 * pi * (x(xa < x & x < xb) - xa) / (xb - xa)).^2;
h(xb < x & x <= xc) = h_recif;
h(xc < x & x < xd) = h_recif - (h_recif - h_ocean) * ...
                     sin(0.5 * pi * (xc - x(xc < x & x < xd)) / (xc - xd)).^2;
disp('Data loaded')

%% Amplitude diff plot
figure
srf1 = surf(x,t(1:end-1),diff(f));
srf1.LineStyle = 'none';

%% Tsunami velocity and amplitude comparison
% Theoretical data
uth = sqrt(g * h(1:end-1));
aAwkb = h.^(1/4);
aBwkb = h.^(-1/4);
aCwkb = h.^(-3/4);

%% Simulations for fixed x
% t = 12000 seconds gives 6 peaks
peak = 3; % choice of the peak to follow
fr = 20; % range of points for interpolation
lim = 0.1; % level above or below 0 for peak detection
[aA,uA,tmaxA] = maxWave2(x,t,f,peak,fr,lim);
[aB,uB,tmaxB] = maxWave2(x2,t2,f2,peak,fr,lim);
[aC,uC,tmaxC] = maxWave2(x3,t3,f3,peak,fr,lim);
aAth = aA(1)/aAwkb(1)*aAwkb;
aBth = aB(1)/aBwkb(1)*aBwkb;
aCth = aC(1)/aCwkb(1)*aCwkb;

%% Fit analysis
figure
hold on
% 456:459 B maximum shift problem
for i = 456:459
    plot(t,f(:,i),'.')
    plot([tmaxB(i) tmaxB(i)],[min(f2(:,i)) max(f2(:,i))],'r')
end
hold off

%% Velocity
% Best result for peak = 1
figure
hold on
plot(x(1:end-1),uth,'k.')
plot(x(1:end-1),uA,'b.')
plot(x(1:end-1),uB,'m.')
plot(x(1:end-1),uC,'r.')
hold off

%% Amplitude
% Better results for peak = 2, 3
figure
hold on
plot(x,aAth,'b')
plot(x,aA,'r')
hold off
figure
hold on
plot(x,aBth,'b')
plot(x,aB,'r')
hold off
figure
hold on
plot(x,aCth,'b')
plot(x,aC,'r')
hold off

%% Live wave simulation A
m = 800; % speed parameter
% theoretical amplitude in black
% numerical amplitude in green
insightWave(x,t,f,aA,aAth,m)
%% B
m = 800;
insightWave(x2,t2,f2,aB,aBth,m)
%% C
m = 800;
insightWave(x3,t3,f3,aC,aCth,m)




%% Tsunami b)
nsimul = 10;
xavar = linspace(xa,xb-20000,nsimul);
%%
for i = 1:nsimul
    cmd = sprintf('wsl %s%s %s %s=%s %s=%g.15 %s=%s', repertoire, executable, input,...
                  'schema', 'B', 'xa', xavar(i), 'output', ['xa=' num2str(xavar(i))]);
    disp(cmd); system(cmd);
end
%%
for i = 1:nsimul
    fichier = ['xa=' num2str(xavar(i))];
    data = load([fichier,'_u.out']);
    x = data(:,1);
    u = data(:,2);
    data = load([fichier,'_E.out']);
    t = data(:,1);
    E = data(:,2);
    data = load([fichier,'_f.out']);
    f = data(:,2:end);
    
    figure('Name',['Analyse de ' fichier])
    subplot(2,2,1)
    plot(x,u)
    grid
    xlabel('x [m]')
    ylabel('u [m/s]')

    subplot(2,2,2)
    plot(t,E)
    grid
    xlabel('t [s]')
    ylabel('E [m^3]')

    subplot(2,2,4)
    pcolor(x,t,f)
    shading interp
    colormap jet
    c = colorbar;
    xlabel('x [m]')
    ylabel('t [s]')
    ylabel(c,'f(x,t) [m]')

    subplot(2,2,3)
    h = plot(x,f(1,:));
    grid
    xlabel('x [m]')
    ylabel('f(x,t) [m]')
    ht = title('t=0 s');
    ylim([min(f(:)),max(f(:))])
end






