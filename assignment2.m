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

%spy(G)
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

%% Question 2 a)
% This problem investigates current flow through an area with a bottleneck by
% using the same finite difference method used in Q1, and applying
% different conductivities to boxes within the area (1 for outside the box,
% 10^-2 for inside). To model the current flow, a voltage was assumed to be
% applied to one side of the area. I decided to set Vo = 10 at x = 0. The
% boxes with lower conductivity are set for 20 < x < 30, and both 0 < y <
% 20, and 30 < y < W, where L and W are the length and width of the area
% respectively and were chosen to both be 50 for simplicity. 
% The following snippet sets up the problem with the constants required.

nw = 50;
nl = 50;
Vo = 10;
sigma_in = 10E-2;
sigma_out = 1;
wanted_solutions = 6;
sigma = zeros(nl, nw);

% build our sigma matrix and plot it
for i = 1:nl
    for j = 1:nw
        % check if inside box
        if (i > 20 && i < 30) && ((j > 1 && j < 20) || (j > 30))
            sigma(i, j) = sigma_in;
        % otherwise, we are outside the box
        else
            sigma(i, j) = sigma_out;
        end
    end
end
surf(sigma)

% initialize and populate the G matrix
G = zeros(nl*nw, nl*nw);
for i = 1:nl
    for j = 1:nw
        n = j + (i - 1) * nw;
        nxm = j + (i - 2) * nw;
        nxp = j + i * nw;
        nyp = (j + 1) + (i - 1) * nw;
        nym = (j - 1) + (i - 1) * nw;
        
        % check if we are at the l = 0 case, apply boundary
        if i == 1
            G(n, n) = Vo;
       
        % we still need to check if we are at an edge to know which
        % terms contribute to G. Corners are most specific cases, we will
        % cover those first, except for x = 0 because that is tied to a
        % boundary condition and will always be Vo
        
        % top right corner, no nxp or nym
        elseif i == nl && j == 1
            G(n, n) = -4*sigma(i, j);
            G(n, nxm) = 1*sigma(i, j);
            G(nyp, n) = 1*sigma(i, j);
        % bottom right corner, no nxp or nyp
        elseif i == nl && j == nw
            G(n, n) = -4*sigma(i, j);
            G(n, nxm) = 1*sigma(i, j);
            G(nym, n) = 1*sigma(i, j);     
        
        % now we cover the cases along the edge. We already know we aren't
        % at a corner at this point
        
        % along right edge, so no nxp term
        elseif i == nl
            G(n, n) = -4*sigma(i, j);
            G(n, nxm) = 1*sigma(i, j);
            G(nyp, n) = 1*sigma(i, j);
            G(nym, n) = 1*sigma(i, j);
        % along top edge, no nym
        elseif j == 1
            G(n, n) = -4*sigma(i, j);
            G(n, nxm) = 1*sigma(i, j);
            G(n, nxp) = 1*sigma(i, j);
            G(nyp, n) = 1*sigma(i, j);
        % along bottom edge, no nyp
        elseif j == 1
            G(n, n) = -4*sigma(i, j);
            G(n, nxm) = 1*sigma(i, j);
            G(n, nxp) = 1*sigma(i, j);
            G(nym, n) = 1*sigma(i, j);
        % otherwise, we're somewhere in the middle, all terms considered
        else      
            G(n, n) = -4*sigma(i, j);
            G(n, nxm) = 1*sigma(i, j);
            G(n, nxp) = 1*sigma(i, j);
            G(nyp, n) = 1*sigma(i, j);
            G(nym, n) = 1*sigma(i, j);  
        end
    end
end

% make sure we have the diagonal correct, and find eigen[vectors/values]
spy(G)
[E, D] = eigs(G, wanted_solutions, 'SM');

% plot and pray, starting with J(x,y)
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
