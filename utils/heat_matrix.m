
model = createpde(); %create starting environment
importGeometry(model, 'BracketTwoHoles.stl'); %set geometry


% steel parameters
lambda = 52 /1e3;          % thermal conductivity [W/(mm·K)]
rho    = 7.8e-6;      % density [kg/mm³]
cp     = 480;         % specific heat capacity [J/(kg·K)]


%internal faces of the two holes
faces = [9,10,11,12,13,14];


h_conv_ext = 750 /1e6;         % convection coefficient [W/(mm²·K)]

%specify heat equation
specifyCoefficients(model, 'm', 0, 'd', rho*cp, 'c', lambda, 'a', 0, 'f', 0);

%apply bounday condition: uniform convection
applyBoundaryCondition(model, 'neumann', 'Face', 1:model.Geometry.NumFaces, 'q', h_conv_ext, 'g', 0);


%Hmax = 2; 
Hmax_faces = Hmax / 2;

%generate mesg, more refined on internal faces of holes
generateMesh(model, ...
    'Hmax', Hmax, ...
    'Hface', {faces, Hmax_faces}, ...
    'Hgrad', 1.2, ...                      
    'GeometricOrder', 'linear');

FEM = assembleFEMatrices(model); 

%assemble A and E. Make E diagonal through lumping.
A = -(FEM.K + FEM.Q);
A = (A + A')/2;
N = size(A, 1);
M_fem = (FEM.M + FEM.M') / 2;
E = spdiags(sum(M_fem, 2),0,N,N); 
L = chol(E, 'lower');

%c represent heat input on internal faces of the two holes
nodes_input = findNodes(model.Mesh, 'region', 'Face', faces);
c = zeros(N, 1);
e_diag = diag(E);                         
c(nodes_input) = e_diag(nodes_input);
c = c / norm(c);


mult = @(v) -L\(A*((L')\v));
clear M_fem;