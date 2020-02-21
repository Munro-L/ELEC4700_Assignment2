%% Question 1

% setup
clearvars
clearvars -GLOBAL
close all
nw = 40;
nl = 60;
Vo = 100;
wanted_solutions = 2;

%% Question 1 a)
% For 1 a), I will be focusing on the 1D case. So width will not be
% considered. Because of this, the G matrix will only have dimensions of nl
% by nl. This also changes the way the G matrix is populated compared to
% the EIPA. Since there is only one dimesnion being considered, the
% centermost n term will only have a -2. In addition to this nym and nyp
% also don't have to be considered. Note, setting on eboundary to zero
% causes problems with the singularity of the matrix. So to get around
% this, the boundary condition at x = L is instead set to something
% comparatively small (like 1) next to Vo (big, 100 in this case). 

G = zeros(nl, nl);
for i = 1:nl
    n = i;
    nxm = (i - 1);
    nxp = (i + 1);
    if i == 1
        G(n, n) = Vo;
    elseif i == nl
        G(n, n) = 1;
    else
        G(n, n) = -2;
        G(n, nxm) = 1;
        G(n, nxp) = 1;
    end
end

spy(G)
[E, D] = eigs(G, wanted_solutions, 'SM');

for solution = 1:wanted_solutions
    map = zeros(nl,1);
    for count = 1:nl
        map(count) = E(count, solution);
    end
    figure(solution)
    plot(1:nl, map)
end

%% Question 1 b)
% The following snippet solves for the case where V = Vo (set to 100) at x
% - 0, and V = 0 at x = L. Note, actually setting a boundary condition to
% zero causes problems for the eigs function like in 1 a).

G = zeros(nl*nw, nl*nw);
for i = 1:nl
    for j = 1:nw
        n = j + (i - 1)*nw;
        nxm = j + (i - 2)*nw;
        nxp = j + i*nw;
        nym = (j - 1) + (i - 1)*nw;
        nyp = (j + 1) + (i - 1)*nw;
        
        % check if we are at a boundary
        if i == 1
            G(n, n) = Vo;
        elseif i == nl
            G(n, n) = 1;        % setting to zero breaks eigs, so we get as close as we can
        elseif j == 1 || j == nw
            G(n, n) = 1;
        else
            % populate G matrix
            G(n, n) = -4;
            G(n, nxm) = 1;
            G(n, nxp) = 1;
            G(nyp, n) = 1;
            G(nym, n) = 1;
        end
    end
end

spy(G)
[E, D] = eigs(G, wanted_solutions, 'SM');

for solution = 1:wanted_solutions
    i = 1;
    j = 1;
    map = zeros(nw,nl);
    for count = 1:nw*nl
        map(i, j) = E(count, solution);
        if j == nw
            j = 1;
            i = i + 1;
        else
            j = j + 1;
        end
    end
    figure(solution)
    surf(map)
    axis([0 nw 0 nl])
end
