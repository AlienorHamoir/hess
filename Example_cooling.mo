within ;
model Example_cooling
  Buildings.Fluid.Movers.SpeedControlled_y pump(redeclare package Medium =
        Buildings.Media.Water, redeclare
      Buildings.Fluid.Movers.Data.Pumps.Wilo.Stratos25slash1to4 per)
    annotation (Placement(transformation(extent={{40,40},{60,60}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow losses
    annotation (Placement(transformation(extent={{-68,-90},{-48,-70}})));
  Modelica.Blocks.Sources.Sine load(
    amplitude=5000,
    f(displayUnit="Hz") = 0.001,
    offset=2000)
    annotation (Placement(transformation(extent={{-100,-90},{-80,-70}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor fuelcell(C=100)
    annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));
  Buildings.Fluid.Sources.Boundary_pT airtowater(
    redeclare package Medium = Buildings.Media.Water,
    use_T_in=true,
    nPorts=1) annotation (Placement(transformation(extent={{-20,40},{0,60}})));
  Buildings.Fluid.FixedResistances.PressureDrop pipe(
    redeclare package Medium = Buildings.Media.Water,
    m_flow_nominal=0.5,
    dp_nominal(displayUnit="bar") = 500000)
    annotation (Placement(transformation(extent={{10,40},{30,60}})));
  Buildings.Fluid.MixingVolumes.MixingVolume heatexchanger(
    redeclare package Medium = Buildings.Media.Water,
    m_flow_nominal=0.5,
    V=0.5,
    nPorts=2) annotation (Placement(transformation(extent={{60,20},{80,0}})));
  Buildings.Fluid.Sources.Boundary_pT sink(redeclare package Medium =
        Buildings.Media.Water, nPorts=1)
    annotation (Placement(transformation(extent={{100,40},{80,60}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor fuelcellTemperature
    annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));
  Buildings.Controls.Continuous.LimPID controller(k=0.1, reverseActing=false)
    annotation (Placement(transformation(extent={{-60,70},{-40,90}})));
  Modelica.Blocks.Sources.Constant setpoint(k=50 + 273.15)
    annotation (Placement(transformation(extent={{-100,70},{-80,90}})));
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(filNam=
        Modelica.Utilities.Files.loadResource(
        "modelica://Buildings/Resources/weatherdata/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.mos"))
    "Weather data reader"
    annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
  Buildings.BoundaryConditions.WeatherData.Bus
      weaBus "Weather data bus" annotation (Placement(transformation(extent={{-60,40},
            {-40,60}}),          iconTransformation(extent={{190,-10},{210,10}})));
equation
  connect(load.y, losses.Q_flow)
    annotation (Line(points={{-79,-80},{-68,-80}}, color={0,0,127}));
  connect(fuelcell.port, losses.port)
    annotation (Line(points={{-40,-60},{-40,-80},{-48,-80}}, color={191,0,0}));
  connect(airtowater.ports[1], pipe.port_a)
    annotation (Line(points={{0,50},{10,50}}, color={0,127,255}));
  connect(pipe.port_b, pump.port_a)
    annotation (Line(points={{30,50},{40,50}}, color={0,127,255}));
  connect(heatexchanger.heatPort, losses.port) annotation (Line(points={{60,10},
          {-20,10},{-20,-80},{-48,-80}}, color={191,0,0}));
  connect(fuelcellTemperature.port, fuelcell.port)
    annotation (Line(points={{-80,-30},{-80,-60},{-40,-60}}, color={191,0,0}));
  connect(controller.y, pump.y)
    annotation (Line(points={{-39,80},{50,80},{50,62}}, color={0,0,127}));
  connect(fuelcellTemperature.T, controller.u_m)
    annotation (Line(points={{-59,-30},{-50,-30},{-50,68}}, color={0,0,127}));
  connect(setpoint.y, controller.u_s)
    annotation (Line(points={{-79,80},{-62,80}}, color={0,0,127}));
  connect(heatexchanger.ports[1], pump.port_b)
    annotation (Line(points={{69,20},{69,50},{60,50}}, color={0,127,255}));
  connect(weaDat.weaBus, weaBus) annotation (Line(
      points={{-60,50},{-50,50}},
      color={255,204,51},
      thickness=0.5), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}},
      horizontalAlignment=TextAlignment.Left));
  connect(airtowater.T_in, weaBus.TDryBul) annotation (Line(points={{-22,54},{
          -34,54},{-34,50.05},{-49.95,50.05}}, color={0,0,127}), Text(
      string="%second",
      index=1,
      extent={{-6,3},{-6,3}},
      horizontalAlignment=TextAlignment.Right));
  connect(sink.ports[1], heatexchanger.ports[2])
    annotation (Line(points={{80,50},{71,50},{71,20}}, color={0,127,255}));
  annotation (
    uses(Buildings(version="11.0.0"), Modelica(version="4.0.0")),
    experiment(
      StartTime=12960000,
      StopTime=13046400,
      __Dymola_Algorithm="Dassl"),
    Diagram(graphics={Rectangle(
          extent={{-100,-16},{-10,-100}},
          lineColor={255,255,255},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid), Text(
          extent={{-44,-18},{-22,-34}},
          textColor={0,0,0},
          textStyle={TextStyle.Bold},
          textString="Fuelcell
HeatPort")}));
end Example_cooling;
