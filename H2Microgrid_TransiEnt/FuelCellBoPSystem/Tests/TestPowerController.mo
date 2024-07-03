within H2Microgrid_TransiEnt.FuelCellBoPSystem.Tests;
model TestPowerController "Testing of power controller used in FC system"

    extends TransiEnt.Basics.Icons.Checkmodel;

  Modelica.Blocks.Sources.Constant Voltage(k=0) annotation (Placement(transformation(extent={{-86,-60},{-66,-40}})));
  Modelica.Blocks.Sources.Constant Power(k=500) annotation (Placement(transformation(extent={{-82,42},{-62,62}})));
  FuelCell.Controller.PowerController powerController annotation (Placement(transformation(extent={{26,2},{46,22}})));
  Modelica.Thermal.FluidHeatFlow.Examples.Utilities.DoubleRamp Load(
    startTime=10,
    interval=10,
    duration_1=100,
    duration_2=100,
    offset=0,
    height_1=700,
    height_2=-400)
                 annotation (Placement(transformation(extent={{-84,8},{-64,28}})));
  Modelica.Thermal.FluidHeatFlow.Examples.Utilities.DoubleRamp Load1(
    startTime=10,
    interval=10,
    duration_1=100,
    duration_2=100,
    offset=0,
    height_1=40,
    height_2=-40)
                 annotation (Placement(transformation(extent={{-84,-26},{-64,-6}})));
  Modelica.Blocks.Math.Division division annotation (Placement(transformation(extent={{-20,-12},{0,8}})));
equation
  connect(Load.y, powerController.P) annotation (Line(points={{-63,18},{27,18}},              color={0,0,127}));
  connect(Load.y, division.u1) annotation (Line(points={{-63,18},{-30,18},{-30,4},{-22,4}},          color={0,0,127}));
  connect(powerController.y, division.u2) annotation (Line(points={{47,12},{52,12},{52,-16},{-30,-16},{-30,-8},{-22,-8}},
                                                                                                                color={0,0,127}));
  connect(division.y, powerController.V_stack) annotation (Line(points={{1,-2},{20,-2},{20,6.6},{27,6.6}},      color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=10,
      Interval=1,
      __Dymola_Algorithm="Dassl"),
    Documentation(info="<html>
<p><span style=\"font-family: Arial;\">Test of advanced power controller, based on input stack voltage and power setpoint</span></p>
</html>"));
end TestPowerController;
