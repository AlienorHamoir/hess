within H2Microgrid_TransiEnt.HybridMicrogrid;
model TestMicrogrid "Testing of hybrid microgrid"

    extends TransiEnt.Basics.Icons.Checkmodel;
      parameter String weather_file = Modelica.Utilities.Files.loadResource("modelica://H2Microgrid_TransiEnt/Resources/weather/USA_CA_Los.Angeles.Intl.AP.722950_TMY3.mos") "Path to weather file";

  H2Microgrid_HP Microgrid_HP annotation (Placement(transformation(extent={{-22,-10},{34,38}})));
  Modelica.Blocks.Sources.Sine BESScommand(
    amplitude=5000,
    f(displayUnit="Hz") = 0.01,
    offset=2000) annotation (Placement(transformation(extent={{-74,28},{-54,48}})));
  Modelica.Blocks.Sources.Sine HESScommand(
    amplitude=3000,
    f(displayUnit="Hz") = 0.0001,
    offset=1000) annotation (Placement(transformation(extent={{-80,-18},{-60,2}})));
equation
  connect(BESScommand.y, Microgrid_HP.P_set_battery) annotation (Line(points={{-53,38},{-32,38},{-32,26.48},{-23.12,26.48}}, color={0,0,127}));
  connect(HESScommand.y, Microgrid_HP.P_set_HESS) annotation (Line(points={{-59,-8},{-32,-8},{-32,4.4},{-23.68,4.4}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=28000,
      Interval=1,
      __Dymola_Algorithm="Dassl"),
    Documentation(info="<html>
<p>Testig before integratio as an FMU and interfacing through Python</p>
</html>"));
end TestMicrogrid;
