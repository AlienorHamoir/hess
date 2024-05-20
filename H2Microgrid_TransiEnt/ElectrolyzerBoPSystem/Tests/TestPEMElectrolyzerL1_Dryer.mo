within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Tests;
model TestPEMElectrolyzerL1_Dryer

  extends TransiEnt.Basics.Icons.Checkmodel;
  import TransiEnt;
  import      Modelica.Units.SI;

  parameter TILMedia.VLEFluidTypes.TILMedia_SplineWater water;
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_gas=simCenter.gasModel4;

public
  parameter SI.Power P_el_n=1e6 "Nominal electrical power of the electrolyzer";
  parameter SI.Power P_el_min=0.05*P_el_n "Minimal electrical power of the electrolyzer";
  parameter SI.Power P_el_max=1.0*P_el_n "Maximal electrical power of the electrolyzer";
  parameter SI.Temperature T_out=273.15+15 "Temperature of the produced hydrogen";
  parameter SI.Efficiency eta_n=0.75 "Nominal efficiency of the electrolyzer";
  parameter SI.Pressure p_out=50e5 "Pressure of the produced hydrogen";

  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_pTxi sink_syngas(
    medium=medium_gas,
    variable_p=false,
    p_const=100000) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={70,-70})));
  Electrolyzer.ElectrolyzerL1System_Dryer electrolyzerSystem(medium=medium_gas) annotation (Placement(transformation(extent={{-30,4},{16,46}})));
  TransiEnt.Components.Boundaries.Electrical.ActivePower.Frequency ElectricGrid_0thOrder annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-68,26})));
  Modelica.Blocks.Sources.Ramp PowerRamp(
    offset=2,
    startTime=0,
    duration=3,
    height=40000) annotation (Placement(transformation(extent={{-88,52},{-68,72}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y=5) annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=180,
        origin={46,42})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{-92,-62},{-72,-42}})));
  inner TransiEnt.ModelStatistics modelStatistics annotation (Placement(transformation(extent={{-58,-62},{-38,-42}})));
  Modelica.Blocks.Sources.Ramp PressureRamp(
    offset=0,
    startTime=0,
    duration=0,
    height=100000) annotation (Placement(transformation(extent={{-56,74},{-36,94}})));
  H2DryingSystem.H2DryingSubsystem h2DryingSubsystem annotation (Placement(transformation(extent={{6,-54},{26,-34}})));
equation
  connect(ElectricGrid_0thOrder.epp,electrolyzerSystem. epp) annotation (Line(
      points={{-58,26},{-38,26},{-38,25},{-30,25}},
      color={0,135,135},
      thickness=0.5));
  connect(PowerRamp.y,electrolyzerSystem. P_el_set) annotation (Line(points={{-67,62},{-7,62},{-7,46.84}}, color={0,0,127}));
  connect(realExpression.y,electrolyzerSystem. m_flow_feedIn) annotation (Line(points={{35,42},{24,42},{24,41.8},{16,41.8}},      color={0,0,127}));
  connect(PressureRamp.y, electrolyzerSystem.pressureIn) annotation (Line(points={{-35,84},{6.8,84},{6.8,48.1}}, color={0,0,127}));
  connect(electrolyzerSystem.gasPortOut, h2DryingSubsystem.gasPortIn) annotation (Line(
      points={{-7,4.21},{-8,4.21},{-8,-44},{6.2,-44}},
      color={255,255,0},
      thickness=1.5));
  connect(h2DryingSubsystem.gasPortOut, sink_syngas.gasPort) annotation (Line(
      points={{26,-44},{54,-44},{54,-70},{60,-70}},
      color={255,255,0},
      thickness=1.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end TestPEMElectrolyzerL1_Dryer;
