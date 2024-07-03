within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.CoolingSystem.HeatPortCooling;
model CoolingModel_Valve

  Buildings.Fluid.Sources.Boundary_pT watertowater(
    redeclare package Medium = Buildings.Media.Water,
    use_p_in=false,
    p=400000,
    use_T_in=false,
    T=292.15,
    nPorts=1) annotation (Placement(transformation(extent={{-62,-10},{-42,10}})));
  Buildings.Fluid.MixingVolumes.MixingVolume heatexchanger(
    redeclare package Medium = Buildings.Media.Water,
    m_flow_nominal=5e-2,
    V=0.002,
    nPorts=2) annotation (Placement(transformation(extent={{46,-54},{66,-74}})));
  Buildings.Fluid.Sources.Boundary_pT sink(redeclare package Medium = Buildings.Media.Water, nPorts=1)
    annotation (Placement(transformation(extent={{88,-10},{68,10}})));
  Buildings.Controls.Continuous.LimPID controller(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=22,
    Ti=10,
    Td=3,
    yMax=1,
    yMin=0,
    Ni=0.95,
    Nd=4.8,
    reverseActing=false)
    annotation (Placement(transformation(extent={{-78,56},{-58,76}})));
  Modelica.Blocks.Interfaces.RealInput T_op annotation (Placement(transformation(extent={{-132,46},{-92,86}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor fuelcellTemperature
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-68,-42})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPortCooling "Heat port for heat exchange with the control volume" annotation (Placement(transformation(extent={{-110,-102},{-90,-82}}), iconTransformation(extent={{-110,-102},{-90,-82}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor cellBuffer(
    C=100,
    T(start=296.15, fixed=true),
    der_T(start=0))             annotation (Placement(transformation(extent={{-28,-56},{-8,-36}})));
  Buildings.Fluid.Actuators.Valves.TwoWayLinear coolingValve(
    redeclare package Medium = Buildings.Media.Water,
    m_flow_nominal=10,
    dpValve_nominal=0.001) annotation (Placement(transformation(extent={{14,-10},{34,10}})));
equation
  connect(sink.ports[1], heatexchanger.ports[1]) annotation (Line(points={{68,0},{55,0},{55,-54}},   color={0,127,255}));
  connect(T_op, controller.u_s) annotation (Line(points={{-112,66},{-80,66}}, color={0,0,127}));
  connect(heatPortCooling, cellBuffer.port) annotation (Line(points={{-100,-92},{-18,-92},{-18,-56}}, color={191,0,0}));
  connect(cellBuffer.port, heatexchanger.heatPort) annotation (Line(points={{-18,-56},{-18,-64},{46,-64}}, color={191,0,0}));
  connect(fuelcellTemperature.T, controller.u_m) annotation (Line(points={{-68,-31},{-68,54}},                                                 color={0,0,127}));
  connect(fuelcellTemperature.port, cellBuffer.port) annotation (Line(points={{-68,-52},{-68,-64},{-18,-64},{-18,-56}},           color={191,0,0}));
  connect(controller.y, coolingValve.y) annotation (Line(points={{-57,66},{24,66},{24,12}}, color={0,0,127}));
  connect(coolingValve.port_b, heatexchanger.ports[2]) annotation (Line(points={{34,0},{50,0},{50,-54},{57,-54}}, color={0,127,255}));
  connect(watertowater.ports[1], coolingValve.port_a) annotation (Line(points={{-42,0},{14,0}}, color={0,127,255}));
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
      __Dymola_Algorithm="Dassl"),
    Documentation(info="<html>
<p><span style=\"font-family: Arial;\">Cooling system model: based on stack operating temperature and cooling flux from the cell, a PID actuates the pump speed, which consumes power.</span></p>
<p><span style=\"font-family: Arial;\">The input temperature of the water boundary depends on the outside environment temperature, to account for the cooling of the external water circuit through a cooling tower with its environment. We assume ideal HEX which allows us to lump the different cooling circuits.</span></p>
<p><br><span style=\"font-family: Arial;\">Model elements are chosen in Modelica and Buildings libraries.</span></p>
<p><span style=\"font-family: Arial;\">Model parameters are adapted for a 5kW Giner electrolyzer, based on [1] Z. Abdin, E. MacA. Gray, and C.J. Webb. Modelling and simulation of a proton exchange membrane (PEM) electrolyzer cell. International Journal of Hydrogen Energy, 40(39):13243-13257, 2015. doi:10.1016/j.ijhydene.2015.07.129.</span></p>
<p><br><br>Tested in the check models &quot;H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.CoolingSubsystem.Test_CoolingModel&quot;</p>
<p><br>Model created by Ali&eacute;nor Hamoir in May 2024.</p>
</html>"));
end CoolingModel_Valve;
