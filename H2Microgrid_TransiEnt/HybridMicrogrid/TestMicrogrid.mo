within H2Microgrid_TransiEnt.HybridMicrogrid;
model TestMicrogrid "Testing of hybrid microgrid"

    extends TransiEnt.Basics.Icons.Checkmodel;

  H2Microgrid_HP Microgrid_HP(
    P_set_FC(start=0),
    SOC_start_battery=0.2,
    P_set_battery(start=0)) annotation (Placement(transformation(extent={{-22,-10},{34,38}})));
  Modelica.Blocks.Sources.Sine BESScommand(
    amplitude=5000,
    f(displayUnit="Hz") = 0.01,
    offset=2000) annotation (Placement(transformation(extent={{-82,64},{-62,84}})));
  Modelica.Blocks.Sources.Sine HESScommand(
    amplitude=3000,
    f(displayUnit="Hz") = 0.0001,
    offset=1000) annotation (Placement(transformation(extent={{-80,6},{-60,26}})));
  Modelica.Blocks.Sources.Pulse pulseFC(
    period=400,
    offset=1,
    startTime=10) annotation (Placement(transformation(extent={{-48,-66},{-28,-46}})));
  Modelica.Blocks.Sources.Pulse pulseEL(
    period=400,
    offset=1,
    startTime=40) annotation (Placement(transformation(extent={{-14,-90},{6,-70}})));
  Modelica.Blocks.Sources.Sine HESScommand1(
    amplitude=3000,
    f(displayUnit="Hz") = 0.0001,
    offset=700)  annotation (Placement(transformation(extent={{-78,-28},{-58,-8}})));
equation
  connect(HESScommand.y, Microgrid_HP.P_set_FC) annotation (Line(points={{-59,16},{-32,16},{-32,6.8},{-23.68,6.8}}, color={0,0,127}));
  connect(HESScommand.y, Microgrid_HP.P_set_battery) annotation (Line(points={{-59,16},{-32,16},{-32,26.48},{-23.12,26.48}}, color={0,0,127}));
  connect(pulseFC.y, Microgrid_HP.state_FC) annotation (Line(points={{-27,-56},{-8,-56},{-8,-10.48}}, color={0,0,127}));
  connect(pulseEL.y, Microgrid_HP.state_EL) annotation (Line(points={{7,-80},{19.44,-80},{19.44,-10.96}}, color={0,0,127}));
  connect(HESScommand1.y, Microgrid_HP.P_set_EL) annotation (Line(points={{-57,-18},{-34,-18},{-34,-2.8},{-23.68,-2.8}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=259200,
      __Dymola_NumberOfIntervals=145,
      __Dymola_Algorithm="Dassl"),
    Documentation(info="<html>
<p>Testig before integratio as an FMU and interfacing through Python</p>
</html>"));
end TestMicrogrid;
