%addpath('/home/theo/opt/openEMS/share/CSXCAD/matlab')
%addpath('/home/theo/opt/openEMS/share/openEMS/matlab')

close all
clear
clc

%PARAMETERS
physical_constants;
unit = 1e-3; % specify everything in mm
MSL_length = 17;
MSL_width = 3.2; % Not 50 Ohms for e_r=4.2 and d = 1.58mm;
substrate_thickness = 1.27;
substrate_epr = 10.8;

f_max = 7e9;

%Filter parameters
w_50 = 1.1;
slot_1 = 0.057;
slot_2 = 0.801;

l_1 = 8.148;
l_2 = 8.399;

port_slot = slot_1+l_1+slot_2+l_2+slot_2+l_1+slot_1;


%SETUP THE FDTD PARAMETERS
FDTD = InitFDTD('NrTS', 1e6);
FDTD = SetGaussExcite( FDTD, f_max/2, f_max/2 );
BC   = {'PML_8' 'PML_8' 'MUR' 'MUR' 'PEC' 'MUR'};
FDTD = SetBoundaryCond( FDTD, BC );

%mesh1.x = SmoothMeshLines( [-port_slot/2 l1-port_slot/2+[resolution/3 -resolution/3*2]/4], resolution/4, 1.5 ,0 );
%mesh2.x = SmoothMeshLines( [(l1+l2)-port_slot/2 (l1+l2+l3)-port_slot/2+[resolution/3 -resolution/3*2]/4 ], resolution/4, 1.5 ,0 );
%mesh3.x = SmoothMeshLines( [(l1+l2+l3+l4)-port_slot/2 (l1+l2+l3+l4+l5)-port_slot/2+[resolution/3 -resolution/3*2]/4 ], resolution/4, 1.5 ,0 );
%mesh.x = SmoothMeshLines( [-MSL_length mesh1.x mesh2.x mesh3.x MSL_length], resolution, 1.5 ,0 );
%SETUP THE MESH
CSX = InitCSX();
resolution = c0/(f_max*sqrt(substrate_epr))/unit /50; % resolution of lambda/50
mesh1.x = SmoothMeshLines( [-port_slot/2 slot_1 - port_slot/2+[resolution/3 -resolution/3*2]/4], resolution, 1.5 ,0 );
mesh2.x = SmoothMeshLines( [(-port_slot/2 +slot_1+l_1+slot_2+l_2+slot_2+l_1)  port_slot/2+[resolution/3 -resolution/3*2]/4], resolution, 1.5 ,0 );
mesh.x = SmoothMeshLines( [-MSL_length mesh1.x mesh2.x MSL_length], resolution, 1.5 ,0 );

%mesh1.y = SmoothMeshLines( [0 MSL_width/2+[-resolution/3 +resolution/3*2]/4], resolution/4 , 1.5 ,0);
mesh.y = SmoothMeshLines( [-1.1*MSL_width  1.1*MSL_width], resolution, 1.5 ,0);
mesh.z = SmoothMeshLines( [linspace(0,substrate_thickness,2) 2*substrate_thickness], 2*resolution);
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
 
portstart = [mesh.x(end), -w_50/2, substrate_thickness];
portstop  = [port_slot/2          ,  w_50/2, 0];
[CSX,port{2}] = AddMSLPort( CSX, 999, 2, 'PEC', portstart, portstop, 0, [0 0 -1], 'MeasPlaneShift',  MSL_length/3 );

##w_50 = 1.1;
##slot_1 = 0.057;
##slot_2 = 0.801;
##
##l_1 = 8.148;
##l_2 = 8.399;
start = [-port_slot/2 + slot_1,  w_50/2, substrate_thickness];
stop  = [ -port_slot/2+l_1,  -w_50/2, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );

start = [-port_slot/2 + slot_1+l_1+slot_2,  w_50/2, substrate_thickness];
stop  = [ -port_slot/2+ slot_1+l_1+slot_2+l_2,  -w_50/2, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop );

start = [-port_slot/2 + slot_1+l_1+slot_2+slot_2+l_2,  w_50/2, substrate_thickness];
stop  = [ -port_slot/2+ slot_1+l_1+slot_2+slot_2+l_2+l_1,  -w_50/2, substrate_thickness];
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
 
plot(f/1e9,20*log10(abs(s11)),'k-','LineWidth',2);
hold on;
grid on;
grid minor;
plot(f/1e9,20*log10(abs(s21)),'r--','LineWidth',2);
legend('S_{11}','S_{21}');
ylabel('S-Parameter (dB)','FontSize',12);
xlabel('frequency (GHz) \rightarrow','FontSize',12);
ylim([-70 2]);
