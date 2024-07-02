within H2Microgrid_TransiEnt.FuelCellBoPSystem.Tests;
model TestPowerController

    extends TransiEnt.Basics.Icons.Checkmodel;

  Modelica.Blocks.Sources.Constant Voltage(k=0) annotation (Placement(transformation(extent={{-86,-60},{-66,-40}})));
  Modelica.Blocks.Sources.Constant Power(k=500) annotation (Placement(transformation(extent={{-86,46},{-66,66}})));
  FuelCell.Controller.PowerController powerController annotation (Placement(transformation(extent={{22,-8},{42,12}})));
  Modelica.Thermal.FluidHeatFlow.Examples.Utilities.DoubleRamp Load(
    startTime=10,
    interval=10,
    duration_1=100,
    duration_2=100,
    offset=0,
    height_1=700,
    height_2=-400)
                 annotation (Placement(transformation(extent={{-26,46},{-6,66}})));
  Modelica.Thermal.FluidHeatFlow.Examples.Utilities.DoubleRamp Load1(
    startTime=10,
    interval=10,
    duration_1=100,
    duration_2=100,
    offset=0,
    height_1=40,
    height_2=-40)
                 annotation (Placement(transformation(extent={{-84,-26},{-64,-6}})));
  Modelica.Blocks.Math.Division division annotation (Placement(transformation(extent={{-22,-30},{-2,-10}})));
equation
  connect(Load.y, powerController.P) annotation (Line(points={{-5,56},{18,56},{18,8},{23,8}}, color={0,0,127}));
  connect(Load.y, division.u1) annotation (Line(points={{-5,56},{18,56},{18,-2},{-24,-2},{-24,-14}}, color={0,0,127}));
  connect(powerController.y, division.u2) annotation (Line(points={{43,2},{48,2},{48,-32},{-24,-32},{-24,-26}}, color={0,0,127}));
  connect(division.y, powerController.V_stack) annotation (Line(points={{-1,-20},{18,-20},{18,-3.4},{23,-3.4}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=10,
      Interval=1,
      __Dymola_Algorithm="Dassl"));
end TestPowerController;
