% Example of L and V curve for smoothing of wood surface data
% Paul Eilers, 2018

% Get the data
y = load('WOODSURF.PRN');
m = length(y);
x = (1:m)';
w = 0 * y + 1;  % Optional weights

% Prepare for smoothing
E = speye(m);
D = diff(E, 2);

% Do series of lambdas
fits = [];
pens = [];
lambdas = 10 .^ (0:0.1:8);
for lambda = lambdas
    W = spdiags(w, 0, m, m);
    z = (W + lambda * D' * D) \ (w .* y);
    fit = log10(sum(w .* (y - z) .^ 2));
    pen = log10(sum((D * z) .^ 2));
    fits = [fits; fit];
    pens = [pens; pen];
end

% Compute V-curve
dfits = diff(fits);
dpens = diff(pens);
v = sqrt(dfits .^ 2 + dpens .^ 2);
nla = length(lambdas);
sel = 2:nla;
midla = sqrt(lambdas(sel) .* lambdas(sel - 1));
[vmin k] = min(v);

% Do final smooth
lambda = midla(k);
z = (W + lambda * D' * D) \ (w .* y);
fitopt = log10(sum(w .* (y - z) .^ 2));
penopt = log10(sum((D * z) .^ 2));

% Make plots
figure(1)
subplot(1, 2, 1)
plot(fits, pens, '.')
xlabel('log10(Fits)')
ylabel('log10(Penalties)')
title('L-curve')
hold on
plot(fits(k), pens(k), '+r', 'MarkerSize', 20)
subplot(1, 2, 2)
plot(log10(midla), v, '.')
xlabel('log10(lambda)')
ylabel('Distance')
title('V-curve')
hold on
plot(log10(lambda), vmin, '+r', 'MarkerSize', 20)
hold off

figure(2) 
plot(x, y)
hold on
plot(x, z, 'r', 'LineWidth', 2)
hold off
title('Sanded wood surface')



shg


