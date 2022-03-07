clear subL subR hcL hcR;

%%% load subiculum
[nel,w3d] = read_vrml('dorsal_SWexported.wrl');

%%% is already split!
subL.Vertices = w3d(1).pts;
subL.Faces = w3d(1).knx;
subR.Vertices = w3d(2).pts;
subR.Faces = w3d(2).knx;

%%% load CA1
[nel,w3d] = read_vrml('field_swexported.wrl');

hcAll.Vertices = w3d(1).pts;
hcAll.Faces = w3d(1).knx;

% not yet split
right = find(hcAll.Vertices(:,1) < 0);
left = find(hcAll.Vertices(:,1) >= 0);

%%% get left CA1
hcL.Vertices = hcAll.Vertices(left,:);

nV = length(left);
fcount = 1; clear hcL.Faces;
for iV = 1:nV
    goodFaces = find(any((hcAll.Faces == left(iV)),2)); % matching faces
    for iF = 1:length(goodFaces)
        v = hcAll.Faces(goodFaces(iF),:); % old vertices in this face
        for iDim = 3:-1:1
            new_ind = find(left == v(iDim));
            
            if ~isempty(new_ind)
                temp(iDim) = new_ind; % find new index of old vertex
            else % face crosses over
                temp(iDim) = NaN;
            end
        end
        
        if ~any(isnan(temp))
            hcL.Faces(fcount,:) = temp;
            fcount = fcount + 1;
        end
    end
    
end

% remove duplicates
%hcL.Faces = sort(hcL.Faces,2);
%hcL.Faces = unique(hcL.Faces,'rows');

%%% get right CA1
hcR.Vertices = hcAll.Vertices(right,:);

nV = length(right);
fcount = 1; clear hcR.Faces;
for iV = 1:nV
    goodFaces = find(any((hcAll.Faces == right(iV)),2)); % matching faces
    for iF = 1:length(goodFaces)
        v = hcAll.Faces(goodFaces(iF),:); % old vertices in this face
        for iDim = 3:-1:1
            new_ind = find(right == v(iDim));
            
            if ~isempty(new_ind)
                temp(iDim) = new_ind; % find new index of old vertex
            else % face crosses over
                temp(iDim) = NaN;
            end
        end
        
        if ~any(isnan(temp))
            hcR.Faces(fcount,:) = temp;
            fcount = fcount + 1;
        end
    end
    
end

% remove duplicates
%hcR.Faces = sort(hcR.Faces,2);
%hcR.Faces = unique(hcR.Faces,'rows');

%%% load acc
[nel,w3d] = read_vrml('accumbens_swexported.wrl');
%%% is already split!
accL.Vertices = w3d(14).pts;
accL.Faces = w3d(14).knx;
accR.Vertices = w3d(13).pts;
accR.Faces = w3d(13).knx;

%%% load brain
[nel,w3d] = read_vrml('brain_swexported.wrl');
br.Vertices = w3d(15).pts;
br.Faces = w3d(15).knx;


%% one color per half (sub)
subL.FaceVertexCData = [1 0 0];
subR.FaceVertexCData = [0 1 0];

figure
h1 = patch(subL);
material metal
hold on
h2 = patch(subR);
shading flat
axis off

set(gcf,'renderer','opengl');
set(gca,'Color',0*[0.9 0.9 0.9]);
lighting phong; light; drawnow;
camproj('perspective')
daspect([1 1 1]);

% set(h1,'FaceAlpha',0.5,'EdgeAlpha',0.5)
%% one color per half (hc)
hcL.FaceVertexCData = [1 0 0];
hcR.FaceVertexCData = [0 1 0];

figure
h1 = patch(hcL);
material metal
hold on
h2 = patch(hcR);
shading flat
axis off

set(gcf,'renderer','opengl');
set(gca,'Color',0*[0.9 0.9 0.9]);
lighting phong; light; drawnow;
camproj('perspective')
daspect([1 1 1]);

% set(h1,'FaceAlpha',0.5,'EdgeAlpha',0.5)

%% try a gradient
dorsalPole = [-5 0 -5]; % for HC
anteriorPole = [-5 0 -5]; % for Nacc
nColorSteps = 100;
c1 = [0 0 1]; c2 = [1 0 0];

% first, compute distances
nPoints = length(subL.Vertices); clear sub_d;
for iPoint = nPoints:-1:1
    %d(iPoint) = sqrt(sum((subL.Vertices(iPoint,:)-dorsalPole).^2)); % cartesian
    sub_d(iPoint) = subL.Vertices(iPoint,2); % y
end

nPoints = length(hcL.Vertices); clear hc_d;
for iPoint = nPoints:-1:1
    %d(iPoint) = sqrt(sum((subL.Vertices(iPoint,:)-dorsalPole).^2)); % cartesian
    hc_d(iPoint) = hcL.Vertices(iPoint,2); % y
end

nPoints = length(accL.Vertices); clear acc_d;
for iPoint = nPoints:-1:1
    %d(iPoint) = sqrt(sum((subL.Vertices(iPoint,:)-dorsalPole).^2)); % cartesian
    acc_d(iPoint) = accL.Vertices(iPoint,3); % y
end
    
% map distances onto colorstep
[~,sub_bin] = histc(sub_d,linspace(min(sub_d),max(sub_d),nColorSteps));
[~,hc_bin] = histc(hc_d,linspace(min(hc_d),max(hc_d),nColorSteps));
[~,acc_bin] = histc(acc_d,linspace(min(acc_d),max(acc_d),nColorSteps));

% generate colors for each step
for iC = 1:3
    cmap(:,iC) = linspace(c1(iC),c2(iC),nColorSteps);
end

% update
subL.FaceVertexCData = cmap(sub_bin,:);
subR.FaceVertexCData = [0 0 0];

hcL.FaceVertexCData = cmap(hc_bin,:);
hcR.FaceVertexCData = [0 0 0];

accL.FaceVertexCData = cmap(acc_bin,:);
accR.FaceVertexCData = [0 0 0];

br.FaceVertexCData = [0 0 0];

figure
subL_h = patch(subL);
shading flat
hold on
subR_h = patch(subR);
shading flat

hcL_h = patch(hcL);
shading flat
hold on
hcR_h = patch(hcR);
shading flat

accL_h = patch(accL);
shading flat
hold on
accR_h = patch(accR);
shading flat

br_h = patch(br);
shading flat
axis off

set(gcf,'renderer','opengl');
set(gca,'Color',0*[0.9 0.9 0.9]);
lighting phong; l = light; drawnow;
camproj('perspective')
daspect([1 1 1]);

set(subR_h,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5)
set(hcR_h,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5)
set(accR_h,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5)

set(br_h,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none','FaceAlpha',0.1)

%% try a line
dorsalPole = [-5 0 -5]; % for HC
anteriorPole = [-5 0 -5]; % for Nacc
nColorSteps = 100;
c1 = [0 0 1]; c2 = [1 0 0];

% first, compute distances
nPoints = length(subL.Vertices); clear sub_d;
for iPoint = nPoints:-1:1
    %d(iPoint) = sqrt(sum((subL.Vertices(iPoint,:)-dorsalPole).^2)); % cartesian
    sub_d(iPoint) = subL.Vertices(iPoint,2); % y
end

nPoints = length(hcL.Vertices); clear hc_d;
for iPoint = nPoints:-1:1
    %d(iPoint) = sqrt(sum((subL.Vertices(iPoint,:)-dorsalPole).^2)); % cartesian
    hc_d(iPoint) = hcL.Vertices(iPoint,2); % y
end

% rostrolateral: x = 1.67, y = -7, z = 2.88
% caudomedial: x = 1.2, y = -6.6, z = 0.42
% sanity check: x = -0.97, y = -6.4, z = 0.54 (contralat cm)
rl = [1.67 -7 2.88];
cm = [1.2 -6.6 0.42];
cmtest = [-0.97 -6.4 0.54];

nPoints = length(accL.Vertices); clear acc_d;
for iPoint = nPoints:-1:1
    %d(iPoint) = sqrt(sum((subL.Vertices(iPoint,:)-dorsalPole).^2)); % cartesian
    %acc_d(iPoint) = accL.Vertices(iPoint,3); % y
    
    % first project point to line
    p1 = projectPointOnLine(accL.Vertices(iPoint,:),cmtest,rl);
    % then get distance from cm
    acc_d(iPoint) = sqrt(sum((p1-cmtest).^2));
    
end
    
% map distances onto colorstep
[~,sub_bin] = histc(sub_d,linspace(min(sub_d),max(sub_d),nColorSteps));
[~,hc_bin] = histc(hc_d,linspace(min(hc_d),max(hc_d),nColorSteps));
[~,acc_bin] = histc(acc_d,linspace(min(acc_d),max(acc_d),nColorSteps));

% generate colors for each step
for iC = 1:3
    cmap(:,iC) = linspace(c1(iC),c2(iC),nColorSteps);
end

% update
subL.FaceVertexCData = cmap(sub_bin,:);
subR.FaceVertexCData = [0 0 0];

hcL.FaceVertexCData = cmap(hc_bin,:);
hcR.FaceVertexCData = [0 0 0];

accL.FaceVertexCData = cmap(acc_bin,:);
accR.FaceVertexCData = [0 0 0];

br.FaceVertexCData = [0 0 0];

figure
subL_h = patch(subL);
shading flat
hold on
subR_h = patch(subR);
shading flat

hcL_h = patch(hcL);
shading flat
hold on
hcR_h = patch(hcR);
shading flat

accL_h = patch(accL);
shading flat
hold on
accR_h = patch(accR);
shading flat

%br_h = patch(br);
%shading flat
axis off

set(gcf,'renderer','opengl');
set(gca,'Color',0*[0.9 0.9 0.9]);
lighting phong; l = light; drawnow;
camproj('perspective')
daspect([1 1 1]);

set(subR_h,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5)
set(hcR_h,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5)
set(accR_h,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5)

%set(br_h,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none','FaceAlpha',0.1)

