within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.CoolingSystem.HeatPortCooling;
model Test_CoolingModel
    extends TransiEnt.Basics.Icons.Checkmodel;

  Buildings.Fluid.Movers.SpeedControlled_y pump(redeclare package Medium = Buildings.Media.Water,
    redeclare H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.CoolingSystem.HeatPortCooling.CoolingPumpDESL per,
    inputType=Buildings.Fluid.Types.InputType.Continuous)
    annotation (Placement(transformation(extent={{-2,70},{18,90}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow losses
    annotation (Placement(transformation(extent={{-70,-92},{-50,-72}})));
  Modelica.Blocks.Sources.Sine load(
    amplitude=5000,
    f(displayUnit="Hz") = 0.001,
    offset=2000)
    annotation (Placement(transformation(extent={{-102,-92},{-82,-72}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor fuelcell(C=100)
    annotation (Placement(transformation(extent={{-52,-62},{-32,-42}})));
  Buildings.Fluid.Sources.Boundary_pT watertowater(
    redeclare package Medium = Buildings.Media.Water,
    use_T_in=false,
    T=292.15,
    nPorts=1) annotation (Placement(transformation(extent={{-76,70},{-56,90}})));
  Buildings.Fluid.FixedResistances.PressureDrop pipe(
    redeclare package Medium = Buildings.Media.Water,
    m_flow_nominal=0.5,
    dp_nominal(displayUnit="bar") = 500000)
    annotation (Placement(transformation(extent={{-40,70},{-20,90}})));
  Buildings.Fluid.MixingVolumes.MixingVolume heatexchanger(
    redeclare package Medium = Buildings.Media.Water,
    m_flow_nominal=0.5,
    V=0.5,
    nPorts=2) annotation (Placement(transformation(extent={{28,-26},{48,-46}})));
  Buildings.Fluid.Sources.Boundary_pT sink(redeclare package Medium = Buildings.Media.Water, nPorts=1)
    annotation (Placement(transformation(extent={{96,70},{76,90}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor fuelcellTemperature
    annotation (Placement(transformation(extent={{-82,-42},{-62,-22}})));
  Buildings.Controls.Continuous.LimPID controller(
    k=0.1,
    Ti=0.01,
    reverseActing=false)
    annotation (Placement(transformation(extent={{-52,30},{-32,50}})));
  Modelica.Blocks.Sources.Constant setpoint(k=50 + 273.15)
    annotation (Placement(transformation(extent={{-92,30},{-72,50}})));
equation
  connect(load.y,losses. Q_flow)
    annotation (Line(points={{-81,-82},{-70,-82}}, color={0,0,127}));
  connect(fuelcell.port,losses. port)
    annotation (Line(points={{-42,-62},{-42,-82},{-50,-82}}, color={191,0,0}));
  connect(watertowater.ports[1], pipe.port_a) annotation (Line(points={{-56,80},{-40,80}}, color={0,127,255}));
  connect(fuelcellTemperature.port,fuelcell. port)
    annotation (Line(points={{-82,-32},{-82,-62},{-42,-62}}, color={191,0,0}));
  connect(controller.y,pump. y)
    annotation (Line(points={{-31,40},{-8,40},{-8,102},{8,102},{8,92}},
                                                        color={0,0,127}));
  connect(fuelcellTemperature.T,controller. u_m)
    annotation (Line(points={{-61,-32},{-42,-32},{-42,28}}, color={0,0,127}));
  connect(setpoint.y,controller. u_s)
    annotation (Line(points={{-71,40},{-54,40}}, color={0,0,127}));
  connect(pump.port_b, heatexchanger.ports[1]) annotation (Line(points={{18,80},{37,80},{37,-26}}, color={0,127,255}));
  connect(pipe.port_b, pump.port_a) annotation (Line(points={{-20,80},{-2,80}}, color={0,127,255}));
  connect(sink.ports[1], heatexchanger.ports[2]) annotation (Line(points={{76,80},{39,80},{39,-26}}, color={0,127,255}));
  connect(losses.port, heatexchanger.heatPort) annotation (Line(points={{-50,-82},{22,-82},{22,-36},{28,-36}}, color={191,0,0}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false), graphics={
                      Rectangle(
          extent={{-102,-18},{-12,-102}},
          lineColor={255,255,255},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid), Text(
          extent={{-46,-20},{-24,-36}},
          textColor={0,0,0},
          textStyle={TextStyle.Bold},
          textString="Fuelcell
HeatPort")}),
    experiment(
      StartTime=12960000,
      StopTime=13470000,
      Interval=1,
      __Dymola_Algorithm="Dassl"),
    Documentation(info="<html>
<p>Testing of the pumped cooling system model.</p>
</html>"));
end Test_CoolingModel;
