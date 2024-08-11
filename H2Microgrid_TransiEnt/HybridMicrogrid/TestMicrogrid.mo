within H2Microgrid_TransiEnt.HybridMicrogrid;
model TestMicrogrid "Testing of hybrid microgrid"

    extends TransiEnt.Basics.Icons.Checkmodel;

  H2Microgrid_HP Microgrid_HP(
    P_set_FC(start=0),
    P_set_battery(start=0)) annotation (Placement(transformation(extent={{-22,-10},{34,38}})));
  Modelica.Blocks.Sources.Sine BESScommand(
    amplitude=5000,
    f(displayUnit="Hz") = 0.01,
    offset=2000) annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
  Modelica.Blocks.Sources.Sine HESScommand(
    amplitude=1000,
    f(displayUnit="Hz") = 0.0001,
    offset=1500) annotation (Placement(transformation(extent={{-80,10},{-60,30}})));
  Modelica.Blocks.Sources.Pulse pulseFC(
    amplitude=1,
    period=400,
    offset=1,
    startTime=10) annotation (Placement(transformation(extent={{-46,-90},{-26,-70}})));
  Modelica.Blocks.Sources.Pulse pulseEL(
    period=400,
    offset=1,
    startTime=40) annotation (Placement(transformation(extent={{0,-90},{20,-70}})));
  Modelica.Blocks.Sources.Sine HESScommand1(
    amplitude=3000,
    f(displayUnit="Hz") = 0.0001,
    offset=700)  annotation (Placement(transformation(extent={{-80,-22},{-60,-2}})));
  Modelica.Blocks.Sources.Constant PowerSet(k=0)   annotation (Placement(transformation(extent={{-50,-52},{-30,-32}})));
  Modelica.Blocks.Sources.Constant PowerSet1(k=0.56)
                                                   annotation (Placement(transformation(extent={{-52,68},{-32,88}})));
equation
  connect(PowerSet.y, Microgrid_HP.state_FC) annotation (Line(points={{-29,-42},{-8,-42},{-8,-10.48}}, color={0,0,127}));
  connect(PowerSet.y, Microgrid_HP.state_EL) annotation (Line(points={{-29,-42},{-8,-42},{-8,-18},{19.44,-18},{19.44,-10.96}}, color={0,0,127}));
  connect(PowerSet1.y, Microgrid_HP.P_set_battery) annotation (Line(points={{-31,78},{-26,78},{-26,34},{-32,34},{-32,26.48},{-23.12,26.48}}, color={0,0,127}));
  connect(HESScommand.y, Microgrid_HP.P_set_FC) annotation (Line(points={{-59,20},{-32,20},{-32,6.8},{-23.68,6.8}}, color={0,0,127}));
  connect(HESScommand.y, Microgrid_HP.P_set_EL) annotation (Line(points={{-59,20},{-32,20},{-32,-2.8},{-23.68,-2.8}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=100,
      Interval=1,
      __Dymola_Algorithm="Dassl"),
    Documentation(info="<html>
<p>Testig before integratio as an FMU and interfacing through Python</p>
</html>"));
end TestMicrogrid;
