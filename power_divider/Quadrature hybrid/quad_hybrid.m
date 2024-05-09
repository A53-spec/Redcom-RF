%addpath('/home/theo/opt/openEMS/share/CSXCAD/matlab')
%addpath('/home/theo/opt/openEMS/share/openEMS/matlab')

close all
clear
clc

%PARAMETERS
physical_constants;
unit = 1e-3; % specify everything in mm
MSL_length = 20;

substrate_thickness = 1.58;
substrate_epr = 4.2;

f_min = 6e9;
f_max = 6e9;

%Filter parameters
w_50 = 3.127;
w_35 = 5.317;

l_lambda_4 = 19;

port_slot = l_lambda_4;

%SETUP THE FDTD PARAMETERS
FDTD = InitFDTD();
FDTD = SetGaussExcite( FDTD, f_max/2, f_max/2 );
BC   = {'PML_8' 'PML_8' 'MUR' 'MUR' 'PEC' 'MUR'};
FDTD = SetBoundaryCond( FDTD, BC );


%SETUP THE MESH
CSX = InitCSX();
resolution = c0/(f_max*sqrt(substrate_epr))/unit /50; % resolution of lambda/50
mesh.x = SmoothMeshLines( [-MSL_length MSL_length], resolution, 1.5 ,0 );
mesh.y = SmoothMeshLines( [-5*w_50 5*w_50], resolution, 1.5 ,0);
mesh.z = SmoothMeshLines( [linspace(0,substrate_thickness,2) 3*substrate_thickness], resolution);
CSX = DefineRectGrid( CSX, unit, mesh );

%ADD SUBSTRATE
CSX = AddMaterial( CSX, 'FR4' );
CSX = SetMaterialProperty( CSX, 'FR4', 'Epsilon', substrate_epr );
start = [mesh.x(1),   mesh.y(1),   0];
stop  = [mesh.x(end), mesh.y(end), substrate_thickness];
CSX = AddBox( CSX, 'FR4', 0, start, stop );


%ADD PORTS

CSX = AddMetal( CSX, 'PEC' );
portstart = [ mesh.x(1), w_50/2+port_slot/2, substrate_thickness];
portstop  = [ -port_slot/2,  -w_50/2+port_slot/2, 0];
[CSX,port{1}] = AddMSLPort( CSX, 999, 1, 'PEC', portstart, portstop, 0, [0 0 -1], 'ExcitePort', true, 'FeedShift', 10*resolution, 'MeasPlaneShift',  MSL_length/3);
 
portstart = [mesh.x(end), w_50/2+port_slot/2, substrate_thickness];
portstop  = [port_slot/2,   -w_50/2+port_slot/2, 0];
[CSX,port{2}] = AddMSLPort( CSX, 999, 2, 'PEC', portstart, portstop, 0, [0 0 -1], 'MeasPlaneShift',  MSL_length/3 );

portstart = [mesh.x(end), w_50/2-port_slot/2, substrate_thickness];
portstop  = [port_slot/2 ,   -w_50/2-port_slot/2, 0];
[CSX,port{3}] = AddMSLPort( CSX, 999, 3, 'PEC', portstart, portstop, 0, [0 0 -1], 'MeasPlaneShift',  MSL_length/3 );
 
portstart = [ mesh.x(1), w_50/2-port_slot/2, substrate_thickness];
portstop  = [ -port_slot/2,  -w_50/2-port_slot/2, 0];
[CSX,port{4}] = AddMSLPort( CSX, 999, 4, 'PEC', portstart, portstop, 0, [0 0 -1],  'MeasPlaneShift',  MSL_length/3 );
##   
%ADD POWER DIVIDER
##w_50 = 3.127;
##w_35 = 5.317;
##
##l_lambda_4 = 19;

start = [-port_slot/2,  l_lambda_4/2-w_50/2, substrate_thickness];
stop  = [ -port_slot/2+w_50,  -l_lambda_4/2+w_50/2, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );

start = [port_slot/2,  l_lambda_4/2-w_50/2, substrate_thickness];
stop  = [ port_slot/2-w_50,  -l_lambda_4/2+w_50/2, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );

start = [-port_slot/2,  l_lambda_4/2-w_50/2, substrate_thickness];
stop  = [ l_lambda_4/2,  l_lambda_4/2-w_50/2+w_35, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );

start = [-port_slot/2,  -l_lambda_4/2+w_50/2, substrate_thickness];
stop  = [ l_lambda_4/2,  -l_lambda_4/2+w_50/2-w_35, substrate_thickness];
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
f = linspace( 1e6, f_max, 1601 );
port = calcPort( port, Sim_Path, f, 'RefImpedance', 50);
 
s11 = port{1}.uf.ref./ port{1}.uf.inc;
s21 = port{2}.uf.ref./ port{1}.uf.inc;
s31 = port{3}.uf.ref./ port{1}.uf.inc;
s41 = port{4}.uf.ref./ port{1}.uf.inc;
 
plot(f/1e9,20*log10(abs(s11)),'k-','LineWidth',2);
hold on;
grid on;
grid minor;
set(gca, "linewidth", 1, "fontsize", 22)
plot(f/1e9,20*log10(abs(s21)),'b-','LineWidth',2);
plot(f/1e9,20*log10(abs(s31)),'r--','LineWidth',2);
plot(f/1e9,20*log10(abs(s41)),'g-','LineWidth',2);
legend('S_{11}','S_{21}','S_{31}','S_{41}');
ylabel('S-Parameter (dB)','FontSize',12);
xlabel('frequency (GHz) \rightarrow','FontSize',12);
ylim([-70 2]);
