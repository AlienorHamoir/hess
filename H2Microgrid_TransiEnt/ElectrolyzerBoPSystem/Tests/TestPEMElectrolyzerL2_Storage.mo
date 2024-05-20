within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Tests;
model TestPEMElectrolyzerL2_Storage "Test of PEM Electrolyzer L2 connection to storage"

  extends TransiEnt.Basics.Icons.Checkmodel;
  import TransiEnt;
  import      Modelica.Units.SI;

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;

//protected
//   function plotResult
//
//   constant String resultFileName = "TestPEMElectrolyzer_L1_Dynamics.mat";
//
//   output String resultFile;
//
//   algorithm
//     clearlog();
//     assert(cd(Modelica.Utilities.System.getEnvironmentVariable(TransiEnt.Basics.Types.WORKINGDIR)), "Error changing directory: Working directory must be set as environment variable with name 'workingdir' for this script to work.");
//     resultFile :=TransiEnt.Basics.Functions.fullPathName(Modelica.Utilities.System.getEnvironmentVariable(TransiEnt.Basics.Types.WORKINGDIR) + "/" + resultFileName);
//     removePlots(false);
//     createPlot(id=1, position={0, 0, 1563, 749}, y={"electrolyzer_0thOrder.dynamics.H_flow_H2", "electrolyzer_1stOrder.dynamics.H_flow_H2","electrolyzer_2ndOrder.dynamics.H_flow_H2"},range={0.0,52.0,-100000.0,1100000.0},erase=false,grid=true,filename=resultFile,colors={{28,108,200},{238,46,47},{0,140,72}});
//   end plotResult;

public
  parameter SI.Power P_el_n=1e6 "Nominal electrical power of the electrolyzer";
  parameter SI.Power P_el_min=0.05*P_el_n "Minimal electrical power of the electrolyzer";
  parameter SI.Power P_el_max=1.0*P_el_n "Maximal electrical power of the electrolyzer";
  parameter SI.Temperature T_out=273.15+15 "Temperature of the produced hydrogen";
  parameter SI.Efficiency eta_n=0.75 "Nominal efficiency of the electrolyzer";
  parameter SI.Pressure p_out=50e5 "Pressure of the produced hydrogen";

  TransiEnt.Components.Boundaries.Electrical.ActivePower.Frequency ElectricGrid_0thOrder annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-64,-8})));
  Modelica.Blocks.Sources.Ramp PowerRamp(
    offset=0,
    startTime=10,
    duration=500,
    height=9e4)   annotation (Placement(transformation(extent={{-86,40},{-66,60}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y=5) annotation (Placement(transformation(extent={{28,10},{48,30}})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{-88,-96},{-68,-76}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_pTxi boundary_pTxi(medium=medium) annotation (Placement(transformation(extent={{82,-72},{62,-52}})));
  StorageSystem.H2StorageSystem_Compressed h2StorageSystem(
   includeHeatTransfer=false,
   medium=medium,
  controlCompressorAfterELY(
      p_paramBefore=true,
      p_paramAfter=false,
      p_afterCompParam=20000000)) annotation (Placement(transformation(extent={{12,-70},{32,-50}})));
  inner TransiEnt.ModelStatistics modelStatistics annotation (Placement(transformation(extent={{-54,-96},{-34,-76}})));
  Electrolyzer.ElectrolyzerL2System electrolyzerL2System(
   medium=medium,
    usePowerPort=true,
    t_overload=900,
    m_flow_start=1e-4,
    P_el_min=1e5,
    k=1e11,
    p_out=5000000) annotation (Placement(transformation(extent={{-14,-18},{6,2}})));
equation
  connect(h2StorageSystem.H2PortOut, boundary_pTxi.gasPort) annotation (Line(
      points={{32.4,-62.6},{34,-62},{62,-62}},
      color={255,255,0},
      thickness=1.5));
  connect(electrolyzerL2System.P_el_set, PowerRamp.y) annotation (Line(points={{-4,2.4},{-4,50},{-65,50}},                                           color={0,127,127}));
  connect(electrolyzerL2System.epp, ElectricGrid_0thOrder.epp) annotation (Line(
      points={{-14,-8},{-54,-8}},
      color={0,135,135},
      thickness=0.5));
  connect(electrolyzerL2System.pressureIn, h2StorageSystem.pressureTank) annotation (Line(points={{2,3},{2,34},{60,34},{60,-48},{38,-48},{38,-54},{32.6,-54}},               color={0,0,127}));
  connect(electrolyzerL2System.m_flow_feedIn, realExpression.y) annotation (Line(points={{6,0},{24,0},{24,6},{56,6},{56,20},{49,20}},         color={0,0,127}));
  connect(electrolyzerL2System.gasPortOut, h2StorageSystem.H2PortIn) annotation (Line(
      points={{-4,-17.9},{-4,-62.6},{12,-62.6}},
      color={255,255,0},
      thickness=1.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end TestPEMElectrolyzerL2_Storage;
