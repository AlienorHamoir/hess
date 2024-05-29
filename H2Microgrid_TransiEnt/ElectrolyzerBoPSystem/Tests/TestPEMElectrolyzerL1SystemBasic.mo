within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Tests;
model TestPEMElectrolyzerL1SystemBasic "Test of PEM Electrolyzer L1 system connected to boundary - no storage, no cooling subsystem and no H2 drying"

  extends TransiEnt.Basics.Icons.Checkmodel;
  import TransiEnt;
  import Modelica.Units.SI;

//protected
  function plotResult

  constant String resultFileName = "TestPEMElectrolyzer_L1_Dynamics.mat";

  output String resultFile;

  algorithm
    clearlog();
    assert(cd(Modelica.Utilities.System.getEnvironmentVariable(TransiEnt.Basics.Types.WORKINGDIR)), "Error changing directory: Working directory must be set as environment variable with name 'workingdir' for this script to work.");
    resultFile :=TransiEnt.Basics.Functions.fullPathName(Modelica.Utilities.System.getEnvironmentVariable(TransiEnt.Basics.Types.WORKINGDIR) + "/" + resultFileName);
    removePlots(false);
    createPlot(id=1, position={0, 0, 1563, 749}, y={"electrolyzer_0thOrder.dynamics.H_flow_H2", "electrolyzer_1stOrder.dynamics.H_flow_H2","electrolyzer_2ndOrder.dynamics.H_flow_H2"},range={0.0,52.0,-100000.0,1100000.0},erase=false,grid=true,filename=resultFile,colors={{28,108,200},{238,46,47},{0,140,72}});
  end plotResult;

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
    offset=2,
    startTime=0,
    duration=3,
    height=100000) annotation (Placement(transformation(extent={{-86,40},{-66,60}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y=100) annotation (Placement(transformation(extent={{28,10},{48,30}})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{-88,-96},{-68,-76}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_pTxi boundary_pTxi(medium=simCenter.gasModel4) annotation (Placement(transformation(extent={{82,-72},{62,-52}})));
  inner TransiEnt.ModelStatistics modelStatistics annotation (Placement(transformation(extent={{-54,-96},{-34,-76}})));
  Modelica.Blocks.Sources.Ramp PressureRamp(
    offset=0,
    startTime=0,
    duration=0,
    height=100000) annotation (Placement(transformation(extent={{-34,66},{-14,86}})));
equation
  connect(ElectricGrid_0thOrder.epp, electrolyzerSystem.epp) annotation (Line(
      points={{-54,-8},{-34,-8},{-34,-9},{-26,-9}},
      color={0,135,135},
      thickness=0.5));
  connect(PowerRamp.y, electrolyzerSystem.P_el_set) annotation (Line(points={{-65,50},{-3,50},{-3,12.84}}, color={0,0,127}));
  connect(realExpression.y, electrolyzerSystem.m_flow_feedIn) annotation (Line(points={{49,20},{54,20},{54,34},{20,34},{20,7.8}}, color={0,0,127}));
  connect(electrolyzerSystem.gasPortOut, boundary_pTxi.gasPort) annotation (Line(
      points={{-3,-29.79},{-4,-29.79},{-4,-62},{62,-62}},
      color={255,255,0},
      thickness=1.5));
  connect(PressureRamp.y, electrolyzerSystem.pressureIn) annotation (Line(points={{-13,76},{10.8,76},{10.8,14.1}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end TestPEMElectrolyzerL1SystemBasic;
