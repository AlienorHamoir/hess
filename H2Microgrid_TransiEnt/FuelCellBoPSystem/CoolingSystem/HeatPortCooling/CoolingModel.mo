within H2Microgrid_TransiEnt.FuelCellBoPSystem.CoolingSystem.HeatPortCooling;
model CoolingModel

  Buildings.Fluid.Movers.SpeedControlled_y pump(redeclare package Medium = Buildings.Media.Water, redeclare H2Microgrid_TransiEnt.FuelCellBoPSystem.CoolingSystem.HeatPortCooling.CoolingPumpDESL     per)
    annotation (Placement(transformation(extent={{10,-12},{34,12}})));
  Buildings.Fluid.Sources.Boundary_pT watertowater(
    redeclare package Medium = Buildings.Media.Water,
    use_T_in=false,
    T=292.15,
    nPorts=1) annotation (Placement(transformation(extent={{-62,-10},{-42,10}})));
  Buildings.Fluid.FixedResistances.PressureDrop pipe(
    redeclare package Medium = Buildings.Media.Water,
    m_flow_nominal=2.76e-2,
    dp_nominal(displayUnit="bar") = 400000)
    annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));
  Buildings.Fluid.MixingVolumes.MixingVolume heatexchanger(
    redeclare package Medium = Buildings.Media.Water,
    m_flow_nominal=2.76e-2,
    V=0.0004,
    nPorts=2) annotation (Placement(transformation(extent={{46,-54},{66,-74}})));
  Buildings.Fluid.Sources.Boundary_pT sink(redeclare package Medium = Buildings.Media.Water, nPorts=1)
    annotation (Placement(transformation(extent={{88,-10},{68,10}})));
  Buildings.Controls.Continuous.LimPID controller(
    controllerType=Modelica.Blocks.Types.SimpleController.PID,
    k=10,
    Ti=0.001,
    Td=100,
    yMax=1,
    yMin=0,
    Ni=0.5,
    Nd=1,
    reverseActing=false)
    annotation (Placement(transformation(extent={{-78,56},{-58,76}})));
  Modelica.Blocks.Interfaces.RealInput T_op annotation (Placement(transformation(extent={{-132,46},{-92,86}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor fuelcellTemperature
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-68,-42})));
  Modelica.Blocks.Interfaces.RealOutput P_coolingPump "Electrical power consumed" annotation (Placement(transformation(extent={{94,24},{114,44}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPortCooling "Heat port for heat exchange with the control volume" annotation (Placement(transformation(extent={{-110,-102},{-90,-82}}), iconTransformation(extent={{-110,-102},{-90,-82}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor cellBuffer(
    C=100,
    T(start=298.15, fixed=true),
    der_T(start=0))             annotation (Placement(transformation(extent={{-28,-56},{-8,-36}})));
equation
  connect(watertowater.ports[1], pipe.port_a) annotation (Line(points={{-42,0},{-30,0}},   color={0,127,255}));
  connect(controller.y,pump. y)
    annotation (Line(points={{-57,66},{22,66},{22,14.4}},
                                                        color={0,0,127}));
  connect(pump.port_b, heatexchanger.ports[1]) annotation (Line(points={{34,0},{58,0},{58,-54},{55,-54}},
                                                                                                   color={0,127,255}));
  connect(pipe.port_b, pump.port_a) annotation (Line(points={{-10,0},{10,0}},   color={0,127,255}));
  connect(sink.ports[1], heatexchanger.ports[2]) annotation (Line(points={{68,0},{57,0},{57,-54}},   color={0,127,255}));
  connect(T_op, controller.u_s) annotation (Line(points={{-112,66},{-80,66}}, color={0,0,127}));
  connect(pump.P, P_coolingPump) annotation (Line(points={{35.2,10.8},{64,10.8},{64,34},{104,34}},     color={0,0,127}));
  connect(heatPortCooling, cellBuffer.port) annotation (Line(points={{-100,-92},{-18,-92},{-18,-56}}, color={191,0,0}));
  connect(cellBuffer.port, heatexchanger.heatPort) annotation (Line(points={{-18,-56},{-18,-64},{46,-64}}, color={191,0,0}));
  connect(fuelcellTemperature.T, controller.u_m) annotation (Line(points={{-68,-31},{-68,54}},                                                 color={0,0,127}));
  connect(fuelcellTemperature.port, cellBuffer.port) annotation (Line(points={{-68,-52},{-68,-64},{-18,-64},{-18,-56}},           color={191,0,0}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false), graphics={
                                          Text(
          extent={{-110,-68},{-88,-84}},
          textColor={0,0,0},
          textStyle={TextStyle.Bold},
          textString="Fuelcell/Electrolyzer
HeatPort")}),
    experiment(
      StartTime=12960000,
      StopTime=13470000,
      Interval=1,
      __Dymola_Algorithm="Dassl"));
end CoolingModel;
