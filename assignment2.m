% setup
clearvars
clearvars -GLOBAL
close all
nw = 40;
nl = 60;
Vo = 100;
wanted_solutions = 4;

%% Question 1 a)
% For 1 a), I will be focusing on the 1D case. So width will not be
% considered. Because of this, the G matrix will only have dimensions of nl
% by nl. This also changes the way the G matrix is populated compared to
% the EIPA. Since there is only one dimesnion being considered, the
% centermost n term will only have a -2. In addition to this nym and nyp
% also don't have to be considered. Note, setting one boundary to zero
% causes problems with the singularity of the matrix and breaks the eigs function. So to get around
% this, the boundary condition at x = L is instead set to something
% comparatively small (like 1) next to Vo (big, 100 in this case). 

% initialize and populate the G matrix
G = zeros(nl, nl);
for i = 1:nl
    n = i;
    nxm = (i - 1);
    nxp = (i + 1);
    
    % check for left boundary condition
    if i == 1
        G(n, n) = Vo;
        
    % check for right boundary condition
    elseif i == nl
        G(n, n) = 1;
    
    % not at boundary, populate matrix normally
    else
        G(n, n) = -2;
        G(n, nxm) = 1;
        G(n, nxp) = 1;
    end
end

% make sure we have a "diagonal-ish" G matrix and call eigs
figure(1)
spy(G)
title('Spy of G matrix')
[E, D] = eigs(G, wanted_solutions, 'SM');

% loop through solutions and plot them
for solution = 1:wanted_solutions
    map = zeros(nl,1);
    for count = 1:nl
        map(count) = E(count, solution);
    end
    figure(solution+1)
    plot(1:nl, map)
    title(sprintf('Q1 a) Solution %d', solution))
end

%% Question 1 b)
% The following snippet solves for the case where V = Vo (set to 100) at x
% = 0, and V = 0 at x = L. Note, actually setting a boundary condition to
% zero causes problems for the eigs function like in 1 a). This part is
% pretty similar to a), the only difference being we have two dimesnions
% like the EIPA. So we now have a nested for loop populating the G matrix
% with nyp and nym now being considered.

% initialize and populate the G matrix
G = zeros(nl*nw, nl*nw);
for i = 1:nl
    for j = 1:nw
        n = j + (i - 1)*nw;
        nxm = j + (i - 2)*nw;
        nxp = j + i*nw;
        nym = (j - 1) + (i - 1)*nw;
        nyp = (j + 1) + (i - 1)*nw;
        
        % check if we are at a boundary, setting to zero breaks eigs so we get as close as we can
        if i == 1
            G(n, n) = Vo;
        elseif i == nl
            G(n, n) = 1;
        elseif j == 1 || j == nw
            G(n, n) = 1;
            
        % not at a boundary, populate G matrix
        else
            G(n, n) = -4;
            G(n, nxm) = 1;
            G(n, nxp) = 1;
            G(nyp, n) = 1;
            G(nym, n) = 1;
        end
    end
end

% check if diagonal-ish and call eigs
figure(wanted_solutions + 2)
spy(G)
title('Spy of G matrix')
[E, D] = eigs(G, wanted_solutions, 'SM');

% loop though solutions and plot them (surface this time)
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
    figure(solution + wanted_solutions + 3)
    surf(map)
    title(sprintf('Q1 b) Solution %d', solution))
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
wanted_solutions = 4;
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
title('Plot of Sigma Mapping')

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
figure(2*wanted_solutions + 4)
spy(G)
title('Spy of G matrix')
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
    figure(solution + 2*wanted_solutions + 5)
    surf(map)
    title(sprintf('Q2 a) J(x,y) Solution %d', solution))
    axis([0 nw 0 nl])
end

% plot and pray, the sequel: V(x,y)
% pretty much the same plots, but we multiply by conductivity
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
    map = map .* sigma;
    figure(solution + 3*wanted_solutions + 5)
    surf(map)
    title(sprintf('Q2 a) V(x,y) Solution %d', solution))
    axis([0 nw 0 nl])
end

% plot and pray, round three: E(x,y)
% pretty much the same plots, but we want the derivative of V(x,y)
% for solution = 1:wanted_solutions
%     i = 1;
%     j = 1;
%     map = zeros(nw,nl);
%     for count = 1:nw*nl
%         map(i, j) = E(count, solution);
%         if j == nw
%             j = 1;
%             i = i + 1;
%         else
%             j = j + 1;
%         end
%     end
%     map = map .* sigma;
%     figure(solution + 4*wanted_solutions + 5)
%     surf(map)
%     title(sprintf('Q2 E(x,y) Solution %d', solution))
%     axis([0 nw 0 nl])
% end


%% Question 2 b)
% In this section, we investigate the effect of mesh density. To do this,
% we set W and L to 100, double what they were in part a).
% The bottleneck was also increased in size, so the geometry stays
% constant. Ultimately, this results in the same shape, but with double the
% points considered in both directions.
%
% Unfortunately, I have to copy and paste the whole code from part a) here 
% because Matlab insists that functions must be at the end of the file. If I used the function, the code
% won't be displayed where I want when I publish, so to the TA marking
% this, I'm sorry for the redundant clutter :( I've commented out some of
% the plots that are less essential like the sigma and spy, which won't
% change much from a) to trim a page from the report.

nw = 100;
nl = 100;
Vo = 10;
sigma_in = 10E-2;
sigma_out = 1;
wanted_solutions = 4;
sigma = zeros(nl, nw);

% build our sigma matrix and plot it
for i = 1:nl
    for j = 1:nw
        % check if inside box
        if (i > 40 && i < 60) && ((j > 1 && j < 40) || (j > 60))
            sigma(i, j) = sigma_in;
        % otherwise, we are outside the box
        else
            sigma(i, j) = sigma_out;
        end
    end
end
%surf(sigma)
%title('Plot of Sigma Mapping')

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
%figure(2*wanted_solutions + 4)
%spy(G)
%title('Spy of G matrix')
[E, D] = eigs(G, wanted_solutions, 'SM');

% plot J(x,y)
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
    figure(solution + 2*wanted_solutions + 5)
    surf(map)
    title(sprintf('Q2 b) J(x,y) Solution %d', solution))
    axis([0 nw 0 nl])
end

% plot and pray, the sequel: V(x,y)
% pretty much the same plots, but we multiply by conductivity
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
    map = map .* sigma;
    figure(solution + 3*wanted_solutions + 5)
    surf(map)
    title(sprintf('Q2 b) V(x,y) Solution %d', solution))
    axis([0 nw 0 nl])
end


%% Question 2 c)
% In this section, we investigate the narrowing of the bottleneck.

