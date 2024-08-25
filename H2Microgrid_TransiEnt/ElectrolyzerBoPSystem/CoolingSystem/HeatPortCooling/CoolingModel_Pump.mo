within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.CoolingSystem.HeatPortCooling;
model CoolingModel_Pump

  Buildings.Fluid.Movers.SpeedControlled_y pump(redeclare package Medium = Buildings.Media.Water, redeclare H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.CoolingSystem.HeatPortCooling.CoolingPumpDESL per)
    annotation (Placement(transformation(extent={{22,-12},{46,12}})));

  Buildings.Fluid.Sources.Boundary_pT watertowater(
    redeclare package Medium = Buildings.Media.Water,
    use_T_in=true,
    T=292.15,
    nPorts=1) annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
  Buildings.Fluid.FixedResistances.PressureDrop pipe(
    redeclare package Medium = Buildings.Media.Water,
    m_flow_nominal=5e-2,
    dp_nominal(displayUnit="bar") = 500000)
    annotation (Placement(transformation(extent={{-18,-10},{2,10}})));
  Buildings.Fluid.MixingVolumes.MixingVolume heatexchanger(
    redeclare package Medium = Buildings.Media.Water,
    m_flow_nominal=5e-2,
    V=0.002,
    nPorts=2) annotation (Placement(transformation(extent={{46,-54},{66,-74}})));
  Buildings.Fluid.Sources.Boundary_pT sink(redeclare package Medium = Buildings.Media.Water, nPorts=1)
    annotation (Placement(transformation(extent={{86,-10},{66,10}})));
  Buildings.Controls.Continuous.LimPID controller(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=1000,
    Ti=0.01,
    Td=3,
    yMax=1,
    yMin=0,
    Ni=0.95,
    Nd=4.8,
    reverseActing=false)
    annotation (Placement(transformation(extent={{-60,60},{-40,80}})));
  Modelica.Blocks.Interfaces.RealInput T_op annotation (Placement(transformation(extent={{-120,50},{-80,90}}), iconTransformation(extent={{-120,50},{-80,90}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor fuelcellTemperature
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-68,-32})));
  Modelica.Blocks.Interfaces.RealOutput P_CoolPump "Electrical power consumed" annotation (Placement(transformation(extent={{86,26},{112,52}}), iconTransformation(extent={{86,26},{112,52}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPortCooling "Heat port for heat exchange with the control volume" annotation (Placement(transformation(extent={{-116,-88},{-82,-54}}),  iconTransformation(extent={{-116,-88},{-82,-54}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor cellBuffer(
    C=10,
    T(fixed=true, start=296.65),
    der_T(start=0))             annotation (Placement(transformation(extent={{-28,-56},{-8,-36}})));
  Modelica.Blocks.Interfaces.RealInput T_env "Prescribed boundary temperature from weather file" annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
equation
  connect(watertowater.ports[1], pipe.port_a) annotation (Line(points={{-30,0},{-18,0}},   color={0,127,255}));
  connect(controller.y,pump. y)
    annotation (Line(points={{-39,70},{34,70},{34,14.4}},
                                                        color={0,0,127}));
  connect(pump.port_b, heatexchanger.ports[1]) annotation (Line(points={{46,0},{58,0},{58,-54},{55,-54}},
                                                                                                   color={0,127,255}));
  connect(pipe.port_b, pump.port_a) annotation (Line(points={{2,0},{22,0}},     color={0,127,255}));
  connect(sink.ports[1], heatexchanger.ports[2]) annotation (Line(points={{66,0},{57,0},{57,-54}},   color={0,127,255}));
  connect(T_op, controller.u_s) annotation (Line(points={{-100,70},{-62,70}}, color={0,0,127}));
  connect(pump.P, P_CoolPump) annotation (Line(points={{47.2,10.8},{62,10.8},{62,39},{99,39}}, color={0,0,127}));
  connect(heatPortCooling, cellBuffer.port) annotation (Line(points={{-99,-71},{-18,-71},{-18,-56}},  color={191,0,0}));
  connect(cellBuffer.port, heatexchanger.heatPort) annotation (Line(points={{-18,-56},{-18,-64},{46,-64}}, color={191,0,0}));
  connect(fuelcellTemperature.T, controller.u_m) annotation (Line(points={{-68,-21},{-68,48},{-50,48},{-50,58}},                               color={0,0,127}));
  connect(fuelcellTemperature.port, cellBuffer.port) annotation (Line(points={{-68,-42},{-68,-62},{-18,-62},{-18,-56}},           color={191,0,0}));
  connect(T_env, watertowater.T_in) annotation (Line(points={{-100,0},{-62,0},{-62,4},{-52,4}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false), graphics={
                                          Text(
          extent={{-100,-78},{-78,-94}},
          textColor={0,0,0},
          textStyle={TextStyle.Bold},
          textString="Fuelcell/Electrolyzer
HeatPort",fontSize=8)}),
    experiment(
      StartTime=12960000,
      StopTime=13470000,
      Interval=1,
      __Dymola_Algorithm="Dassl"),
    Documentation(info="<html>
<p>Cooling system model: based on stack operating temperature and cooling flux from the cell, a PID actuates the pump speed, which consumes power. </p>
<p>The input temperature of the water boundary depends on the outside environment temperature, to account for the cooling of the external water circuit through a cooling tower with its environment. We assume ideal HEX which allows us to lump the different cooling circuits.</p>
<p>Model elements are chosen in Modelica and Buildings libraries.</p><p>Model parameters are adapted for a 5kW Giner electrolyzer, based on [1] Z. Abdin, E. MacA. Gray, and C.J. Webb. Modelling and simulation of a proton exchange membrane (PEM) electrolyzer cell. International Journal of Hydrogen Energy, 40(39):13243-13257, 2015. doi:10.1016/j.ijhydene.2015.07.129. </p>
<p><br>Tested in the check models &quot;H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.CoolingSubsystem.Test_CoolingModel&quot;</p>
<p><br>Model created by Ali&eacute;nor Hamoir in May 2024.</p>
</html>"));
end CoolingModel_Pump;
