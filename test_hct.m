close all
clear all

A = load('hct.pts');
x = A(2:end, 2);
y = A(2:end, 3);
tri = delaunay(x, y);
nfig = 0;
nfig = nfig + 1;
figure(nfig)
triplot(tri, x, y);

num = 1:size(tri, 1);
tri = [num', tri];
dim = [size(tri, 1), zeros(1, size(tri, 2) - 1)];
tri = [dim; tri];
dlmwrite('hct.tri', tri);

ntestx = 31;
ntesty = 31;
alpha = 0.0;
beta = 3.0;
gamma = 0.0;
delta = 3.0;
fid = fopen('grille.don', 'w');
fprintf(fid, '%i %i\n', ntestx, ntesty);
fprintf(fid, '%f %f %f %f\n', alpha, beta, gamma, delta);
fclose(fid);

system('./a.out');

res = load('S.res');

x = linspace(alpha, beta, ntestx);
y = linspace(gamma, delta, ntesty);
[X, Y] = meshgrid(x, y);

nfig = nfig + 1;
figure(nfig)
S = reshape(res(:, 1), [ntestx, ntesty])';
surf(X, Y, S)
title('HCT')


nfig = nfig + 1;
figure(nfig)
erreur= reshape(res(:, 2), [ntestx, ntesty])';
surf(X, Y, erreur)
title('Erreur')

