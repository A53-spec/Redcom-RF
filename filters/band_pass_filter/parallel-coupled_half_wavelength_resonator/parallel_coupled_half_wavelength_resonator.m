%addpath('/home/theo/opt/openEMS/share/CSXCAD/matlab')
%addpath('/home/theo/opt/openEMS/share/openEMS/matlab')

close all
clear
clc

%PARAMETERS
physical_constants;
unit = 1e-3; % specify everything in mm
MSL_length = 15;
substrate_thickness = 0.635;
substrate_epr = 10.8;

f_min = 7e9;
f_max = 12e9;

%Filter parameters
w_50 = 0.59;
w1 = 0.385;
w2 = 0.575;
w3 = 0.595;
l1 = 2.852;
l2 = 2.772;
l3 = 2.756;
slot1 = 0.161;
slot2 = 0.540;
slot3 = 0.730;

port_slot = 2*(l1 + l2 + l3);
%filter_width = w_50/2+slot1+w2+slot2+w3+slot3+w3+slot3+w3+slot2+w2+slot1;
filter_width = w_50+slot1+w2+slot2+w3+slot3+w3+slot3+w3+slot2+w2+slot1;

%SETUP THE FDTD PARAMETERS
%FDTD = InitFDTD();
FDTD = InitFDTD('NrTS', 1e9, 'EndCriteria', 1e-7);
FDTD = SetGaussExcite( FDTD, (f_min+f_max)/2, (f_max-f_min)/2 );
BC   = {'PML_8' 'PML_8' 'MUR' 'MUR' 'PEC' 'MUR'};
FDTD = SetBoundaryCond( FDTD, BC );


%SETUP THE MESH
CSX = InitCSX();
resolution = c0/(f_max*sqrt(substrate_epr))/unit /50; % resolution of lambda/50
mesh1.x = SmoothMeshLines( [-port_slot/2-w_50/6 -port_slot/2+w_50/6], resolution/6, 1.5 ,0 );
mesh2.x = SmoothMeshLines( [-port_slot/2-w_50/6+l1 -port_slot/2+w_50/6+l1], resolution/6, 1.5 ,0 );
mesh3.x = SmoothMeshLines( [-port_slot/2-w_50/6+l1+l2 -port_slot/2+w_50/6+l1+l2], resolution/6, 1.5 ,0 );
mesh4.x = SmoothMeshLines( [-port_slot/2-w_50/6+l1+l2+l3 -port_slot/2+w_50/6+l1+l2+l3], resolution/6, 1.5 ,0 );
mesh.x = SmoothMeshLines( [-MSL_length -mesh1.x mesh1.x -mesh2.x mesh2.x -mesh3.x mesh3.x mesh4.x MSL_length], resolution, 1.5 ,0 );
%mesh.x = SmoothMeshLines( [-MSL_length  MSL_length], resolution/2, 1.5 ,0 );

mesh1.y = SmoothMeshLines( [-w_50/2+w1+w_50/2 -w_50/2-slot1-w2], resolution/8, 1.5 ,0 );
mesh2.y = SmoothMeshLines( [-w_50/2-slot1-w2-slot2 -w_50/2-slot1-w2-slot2-w3], resolution/8, 1.5 ,0 );
mesh3.y = SmoothMeshLines( [-w_50/2-slot1-w2-slot2-w3-slot3 -w_50/2-slot1-w2-slot2-w3-slot3-w3], resolution/8, 1.5 ,0 );
mesh4.y = SmoothMeshLines( [-w_50/2-slot1-w2-slot2-w3-slot3-w3-slot3 -w_50/2-slot1-w2-slot2-w3-slot3-w3-slot3-w3], resolution/8, 1.5 ,0 );
mesh5.y = SmoothMeshLines( [-w_50/2-slot1-w2-slot2-w3-slot3-w3-slot3-w3-slot2 -w_50/2-slot1-w2-slot2-w3-slot3-w3-slot3-w3-slot2-w2-slot1-w1-w_50/2], resolution/8, 1.5 ,0 );
mesh.y = SmoothMeshLines( [-(filter_width+2*w_50) mesh1.y mesh2.y mesh3.y mesh4.y mesh5.y 2*w_50], resolution, 1.5 ,0);
mesh.z = SmoothMeshLines( [linspace(0,substrate_thickness,2) 2*substrate_thickness], 4*resolution);
CSX = DefineRectGrid( CSX, unit, mesh );

%ADD SUBSTRATE
CSX = AddMaterial( CSX, 'FR4' );
CSX = SetMaterialProperty( CSX, 'FR4', 'Epsilon', substrate_epr );
start = [mesh.x(1),   mesh.y(1),   0];
stop  = [mesh.x(end), mesh.y(end), substrate_thickness];
CSX = AddBox( CSX, 'FR4', 0, start, stop );


%ADD PORTS

CSX = AddMetal( CSX, 'PEC' );
portstart = [ mesh.x(1), -w_50/2, substrate_thickness];
portstop  = [ -port_slot/2,  w_50/2, 0];
[CSX,port{1}] = AddMSLPort( CSX, 999, 1, 'PEC', portstart, portstop, 0, [0 0 -1], 'ExcitePort', true, 'FeedShift', 10*resolution, 'MeasPlaneShift',  MSL_length/3);
 
portstart = [mesh.x(end), -w_50/2-filter_width, substrate_thickness];
portstop  = [port_slot/2          ,  w_50/2-filter_width, 0];
[CSX,port{2}] = AddMSLPort( CSX, 999, 2, 'PEC', portstart, portstop, 0, [0 0 -1], 'MeasPlaneShift',  MSL_length/3 );
    
%Filter parameters
##w_50 = 0.59;
##w1 = 0.385;
##w2 = 0.575;
##w3 = 0.595;
##l1 = 2.852;
##l2 = 2.772;
##l3 = 2.756;
##slot1 = 0.161;
##slot2 = 0.540;
##slot3 = 0.730;
start = [-port_slot/2,  -w_50/2, substrate_thickness];
stop  = [ -port_slot/2+l1,  -w_50/2+w1, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );
start = [-port_slot/2,  -w_50/2 - slot1-w1, substrate_thickness];
stop  = [ -port_slot/2+l1,  -w_50/2 - slot1, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );

start = [-port_slot/2+l1,  -w_50/2 - slot1, substrate_thickness];
stop  = [ -port_slot/2+l1+l2,  -w_50/2 - slot1-w2, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );
start = [-port_slot/2+l1,  -w_50/2 - slot1-w2-slot2, substrate_thickness];
stop  = [ -port_slot/2+l1+l2,   -w_50/2 - slot1-w2-slot2-w2, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );

start = [-port_slot/2+l1+l2,  -w_50/2 - slot1-w2-slot2, substrate_thickness];
stop  = [ -port_slot/2+l1+l2+l3,   -w_50/2 - slot1-w2-slot2-w3, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );
start = [-port_slot/2+l1+l2,  -w_50/2 - slot1-w2-slot2-w3-slot3, substrate_thickness];
stop  = [ -port_slot/2+l1+l2+l3,   -w_50/2 - slot1-w2-slot2-w3-slot3-w3, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );

start = [-port_slot/2+l1+l2+l3,  -w_50/2 - slot1-w2-slot2-w3-slot3, substrate_thickness];
stop  = [ -port_slot/2+l1+l2+l3+l3,   -w_50/2 - slot1-w2-slot2-w3-slot3-w3, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );
start = [-port_slot/2+l1+l2+l3,  -w_50/2 - slot1-w2-slot2-w3-slot3-w3-slot3, substrate_thickness];
stop  = [ -port_slot/2+l1+l2+l3+l3,   -w_50/2 - slot1-w2-slot2-w3-slot3-w3-slot3-w3, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );

start = [-port_slot/2+l1+l2+l3+l3,  -w_50/2 - slot1-w2-slot2-w3-slot3-w3-slot3-w3, substrate_thickness];
stop  = [ -port_slot/2+l1+l2+l3+l3+l2,   -w_50/2 - slot1-w2-slot2-w3-slot3-w3-slot3-w3+w2, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );
start = [-port_slot/2+l1+l2+l3+l3,  -w_50/2 - slot1-w2-slot2-w3-slot3-w3-slot3-w3-slot2-w2, substrate_thickness];
stop  = [ -port_slot/2+l1+l2+l3+l3+l2,   -w_50/2 - slot1-w2-slot2-w3-slot3-w3-slot3-w3-slot2, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );

start = [-port_slot/2+l1+l2+l3+l3+l2,  -w_50/2 - slot1-w2-slot2-w3-slot3-w3-slot3-w3-slot2-w2, substrate_thickness];
stop  = [ -port_slot/2+l1+l2+l3+l3+l2+l1,   -w_50/2 - slot1-w2-slot2-w3-slot3-w3-slot3-w3-slot2-w2+w1, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );
start = [-port_slot/2+l1+l2+l3+l3+l2,  -w_50/2 - slot1-w2-slot2-w3-slot3-w3-slot3-w3-slot2-w2-slot1-w1, substrate_thickness];
stop  = [ -port_slot/2+l1+l2+l3+l3+l2+l1,   -w_50/2 - slot1-w2-slot2-w3-slot3-w3-slot3-w3-slot2-w2-slot1, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );
%RUN OPENEMS
Sim_Path = 'tmp';
Sim_CSX = 'msl.xml';
 
[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder
 
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );
RunOpenEMS( Sim_Path, Sim_CSX );

%POST PROCESSING
close all
f = linspace( f_min, f_max, 1601 );
port = calcPort( port, Sim_Path, f, 'RefImpedance', 50);
 
s11 = port{1}.uf.ref./ port{1}.uf.inc;
s21 = port{2}.uf.ref./ port{1}.uf.inc;
 
plot(f/1e9,20*log10(abs(s11)),'k-','LineWidth',2);
hold on;
grid on;
grid minor;
plot(f/1e9,20*log10(abs(s21)),'r--','LineWidth',2);
legend('S_{11}','S_{21}');
ylabel('S-Parameter (dB)','FontSize',12);
xlabel('frequency (GHz) \rightarrow','FontSize',12);
ylim([-70 2]);
