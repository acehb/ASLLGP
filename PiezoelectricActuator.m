function [yd,ym]=PiezoelectricActuator(L,H,d31,V,fidelity)
%Inputs L:beam length;H:beam height;d31:piezoelectric d_31 strain coefficient;V:voltage;fidelity can be 1 or 2, fidelity=1:low-fidelity simulator,fidelity=2:high-fidelity simulator.
%Outputs yd:tip deflection of the beam;ym: mass of the beam.
N = 3;
model = createpde(N);
relPermittivity = 12;
H2 = H/2; % height of each layer in meters 
topLayer = [3 4 0 L L 0 0 0 H2 H2];
bottomLayer = [3 4 0 L L 0 -H2 -H2 0 0];
gdm = [topLayer; bottomLayer]';
g = decsg(gdm, 'R1+R2', ['R1'; 'R2']');
% Create a geometry entity and append to the PDE model.
geometryFromEdges(model,g);
% The material in both layers of the beam is Polyvinylidene Fluoride (PVDF), a thermoplastic polymer with piezoelectric behavior.
E = 2.0e9; % Elastic modulus, N/m^2
NU = 0.29; % Poisson's ratio
G = 0.775e9; % Shear modulus, N/m^2
d33 = -3.0e-11;
permittivityFreeSpace = 8.854187817620e-12; % F/m
C11 = E/(1-NU^2);
C12 = NU*C11;
c2d = [C11 C12 0; C12 C11 0; 0 0 G];
pzeD = [0 d31; 0 d33; 0 0];
pzeE = c2d*pzeD;
D_const_stress = [relPermittivity 0; 0 relPermittivity]*permittivityFreeSpace;
% Convert dielectric matrix from constant stress to constant strain
D_const_strain = D_const_stress - pzeD'*pzeE;
c11 = [c2d(1,1) c2d(1,3) c2d(3,3)];
c12 = [c2d(1,3) c2d(1,2); c2d(3,3) c2d(2,3)];
c22 = [c2d(3,3) c2d(2,3) c2d(2,2)];
c13 = [pzeE(1,1) pzeE(1,2); pzeE(3,1) pzeE(3,2)];
c23 = [pzeE(3,1) pzeE(3,2); pzeE(2,1) pzeE(2,2)];
c33 = [D_const_strain(1,1) D_const_strain(2,1) D_const_strain(2,2)];
ctop = [c11(:); c12(:); c22(:); -c13(:); -c23(:); -c33(:)];
cbot = [c11(:); c12(:); c22(:); c13(:); c23(:); -c33(:)];
f = [0 0 0]';
specifyCoefficients(model, 'm', 0,'d', 0,'c', ctop, 'a', 0, 'f', f,'Face',2);
specifyCoefficients(model, 'm', 0,'d', 0,'c', cbot, 'a', 0, 'f', f,'Face',1);
if(fidelity==1) %LF
    hmax = min([L H2])/2; %min([L H2])/k, can change k over¡¾2£¬infinity). must make sure LF is less time-consuming than HF.  
    generateMesh(model,'Hmax',hmax,'GeometricOrder','linear','MesherVersion','R2013a');    
elseif(fidelity==2) %HF
    hmax = min([L H2])/2.5; %min([L H2])/k, make sure k=2.5 and k=2*2.5=5 give nearly same results.
    generateMesh(model,'Hmax',hmax,'GeometricOrder','quadratic','MesherVersion','R2013a');    
end

% Boundary Condition Definition
topedge=[]; bottomedge=[]; leftedges=[];
for i=1:model.Geometry.NumEdges
    I=findNodes(model.Mesh,'region','Edge',i);
    if(all(model.Mesh.Nodes(2,I)>(H2*0.9999)))
        topedge=[topedge; i];
    elseif(all(model.Mesh.Nodes(2,I)<(-H2*0.9999)))
        bottomedge=[bottomedge; i];    
    elseif(all(model.Mesh.Nodes(1,I)<(L*0.0001)))
        leftedges=[leftedges;i];
    end
end
if((length(topedge)~=1)||(length(bottomedge)~=1)||(length(leftedges)~=2))
    disp('error')
    return
end
voltTop = applyBoundaryCondition(model,'mixed','Edge',topedge,'u',V,'EquationIndex',3);
voltBot = applyBoundaryCondition(model,'mixed','Edge',bottomedge,'u',0,'EquationIndex',3);
clampLeft = applyBoundaryCondition(model,'mixed','Edge',leftedges,'u',[0 0],'EquationIndex',1:2);

result = solvepde(model);
rs = result.NodalSolution;
yd= min(rs(:,2));
ym= H*L;

