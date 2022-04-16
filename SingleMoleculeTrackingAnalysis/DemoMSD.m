SPACE_UNITS = 'µm';
TIME_UNITS = 's';

N_PARTICLES = 10;
N_TIME_STEPS = 100;
N_DIM = 2; % 2D

% Typical values taken from studies of proteins diffusing in membranes:
% Diffusion coefficient
D  = 1.0e-3; % µm^2/s
% Time step between acquisition; fast acquisition!
dT = 0.05; % s,

% Area size, just used to disperse particles in 2D. Has no impact on
% analysis.
SIZE = 2; % µm

%Following Einstein equation on the previous page, displacements follow a Gaussian PDF with standard deviation given by:

k = sqrt(N_DIM * D * dT);

%Let's generate the tracks. @msdanalyzer imposes that the tracks you give to it are formatted in the following way: [ Ti Xi Yi ...]. So if we generate a track with 50 measurements in a 2D diffusion problem, we must generate a 50 x 3 double array per particle.

tracks = cell(N_PARTICLES, 1);

for i = 1 : N_PARTICLES

    % Time
    time = (0 : N_TIME_STEPS-1)' * dT;

    % Initial position
    X0 = SIZE .* rand(1, N_DIM);

    % Integrate uncorrelated displacement
    dX = k * randn(N_TIME_STEPS, N_DIM);
    dX(1, :) = X0;
    X = cumsum(dX, 1);

    % Store
    tracks{i} = [time X];

end
clear i X dX time X0


pause;
%To instantiate the analyzer, we must first provide it with the dimensionality of the problem, the space units and the time units. The two later arguments are just used for convenience.

ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);

%To pass the tracks we created above, we must use obj = obj.doSomething notation, for @msdanalyzer is a per-value class. That is: if you do not catch the object returned, your modifications are lost. More precisely:

% This does not work:
ma.addAll(tracks);

% This works:
ma = ma.addAll(tracks);

% Indeed:
disp(ma)

ma.plotTracks;
ma.labelPlotTracks;

ma = ma.computeMSD;
ma.msd


t = (0 : N_TIME_STEPS)' * dT;

%The calculation of all possible delays is the following:

[T1, T2] = meshgrid(t, t);
all_delays = unique( abs(T1 - T2) );

%We expect to have N different delays (0 included), but got:

fprintf('Found %d different delays.\n', numel(all_delays));

disp(all_delays(1:9));

fprintf('For %d time-points, found %d different delays.\n', N_TIME_STEPS, size( ma.msd{1}, 1 ) );

%To plot the resulting individual MSD curves, use:
figure
ma.plotMSD;

cla
ma.plotMeanMSD(gca, true)


mmsd = ma.getMeanMSD;
t = mmsd(:,1);
x = mmsd(:,2);
dx = mmsd(:,3) ./ sqrt(mmsd(:,4));
errorbar(t, x, dx, 'k');


[fo, gof] = ma.fitMeanMSD;
plot(fo)
ma.labelPlotMSD;
legend off


ma = ma.fitMSD;
good_enough_fit = ma.lfit.r2fit > 0.8;
Dmean = mean( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;

fprintf('Estimation of the diffusion coefficient from linear fit of the MSD curves:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
    Dmean, Dstd, sum(good_enough_fit));

pause;

% Retrieve instantaneous velocities, per track
 trackV = ma.getVelocities;

 % Pool track data together
 TV = vertcat( trackV{:} );

 % Velocities are returned in a N x (nDim+1) array: [ T Vx Vy ...]. So the
 % velocity vector in 2D is:
 V = TV(:, 2:3);

 % Compute diffusion coefficient
varV = var(V);
mVarV = mean(varV); % Take the mean of the two estimates
Dest = mVarV / 2 * dT;

fprintf('Estimation from velocities histogram:\n')
fprintf('D = %.3g %s, real value was %.3g %s\n', ...
    Dest, [SPACE_UNITS '²/' TIME_UNITS], D, [SPACE_UNITS '²/' TIME_UNITS]);