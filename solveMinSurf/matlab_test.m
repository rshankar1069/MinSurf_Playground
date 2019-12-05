numberOfPDE = 1;
model = createpde(numberOfPDE);
% 
% pderect([0 1 0 1]);
% g = decsg(R);

geometryFromEdges(model,@squareg);

pdegplot(model,'EdgeLabels','on'); 
axis equal
title 'Geometry with Edge Labels';

a = 0;
f = 0;
cCoef = @(region,state) 1./sqrt(1+state.ux.^2 + state.uy.^2);
specifyCoefficients(model,'m',0,'d',0,'c',cCoef,'a',a,'f',f);

bcMatrix = @(region,~)region.x.^2+region.y.^2;
applyBoundaryCondition(model,'dirichlet',...
                             'Edge',1:model.Geometry.NumEdges,...
                             'u',bcMatrix);
                         
                         generateMesh(model,'Hmax',0.1);
figure; 
pdemesh(model); 
axis equal

model.SolverOptions.ReportStatistics = 'on';
result = solvepde(model);

uu = result.NodalSolution;

figure; 
pdeplot(model,'XYData',uu,'ZData',uu, 'Mesh', 'on');
xlabel 'x'
ylabel 'y'
zlabel 'u(x,y)'
title 'Minimal Surface'