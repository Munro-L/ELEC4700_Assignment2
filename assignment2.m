%% Question 1

% setup
clearvars
clearvars -GLOBAL
close all
nw = 40;
nl = 60;
Vo = 100;
G = zeros(nl*nw, nl*nw);
wanted_solutions = 2;

%% Question 1 a)
% The following snippet solves for the case where V = Vo (set to 100) at x
% - 0, and V = 0 at x = L. Note, actually setting a boundary condition to
% zero causes problems for the eigs function, so instead, V is being set to
% 1 at x = L. This is the reason why Vo was set large.

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