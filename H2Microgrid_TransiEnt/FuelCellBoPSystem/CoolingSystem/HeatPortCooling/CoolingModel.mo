within H2Microgrid_TransiEnt.FuelCellBoPSystem.CoolingSystem.HeatPortCooling;
model CoolingModel


  parameter Real k_p=0.5 "gain, cooling system PID proportional control - 1050 when opposite sign convention with PID";
  parameter Modelica.Units.SI.Time tau_i=0.1 "1/tau_i for cooling system PID integrator gain";
  parameter Real N_i=0.5 "gain of anti-windup compensation ";

  Buildings.Fluid.Movers.SpeedControlled_y pump(redeclare package Medium = Buildings.Media.Water, redeclare Buildings.Fluid.Movers.Data.Pumps.Wilo.TopS40slash10 per)
    annotation (Placement(transformation(extent={{10,-12},{34,12}})));
  Buildings.Fluid.Sources.Boundary_pT watertowater(
    redeclare package Medium = Buildings.Media.Water,
    use_T_in=true,
    T=292.15,
    nPorts=1) annotation (Placement(transformation(extent={{-62,-10},{-42,10}})));
  Buildings.Fluid.FixedResistances.PressureDrop pipe(
    redeclare package Medium = Buildings.Media.Water,
    m_flow_nominal=3e-2,
    dp_nominal(displayUnit="bar") = 100000)
    annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));
  Buildings.Fluid.MixingVolumes.MixingVolume heatexchanger(
    redeclare package Medium = Buildings.Media.Water,
    m_flow_nominal=3e-2,
    V=0.006,
    nPorts=2) annotation (Placement(transformation(extent={{46,-54},{66,-74}})));
  Buildings.Fluid.Sources.Boundary_pT sink(redeclare package Medium = Buildings.Media.Water, nPorts=1)
    annotation (Placement(transformation(extent={{88,-10},{68,10}})));
  Buildings.Controls.Continuous.LimPID controller(
    controllerType=Modelica.Blocks.Types.SimpleController.P,
    k=k_p,
    Ti=tau_i,
    Td=0.01,
    yMax=1,
    yMin=0,
    Ni=N_i,
    Nd=1,
    reverseActing=false)
    annotation (Placement(transformation(extent={{-78,56},{-58,76}})));
  Modelica.Blocks.Interfaces.RealInput T_op annotation (Placement(transformation(extent={{-120,40},{-80,80}}), iconTransformation(extent={{-120,40},{-80,80}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor fuelcellTemperature
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-68,-42})));
  Modelica.Blocks.Interfaces.RealOutput P_coolingPump "Electrical power consumed" annotation (Placement(transformation(extent={{80,20},{120,60}}), iconTransformation(extent={{80,20},{120,60}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPortCooling "Heat port for heat exchange with the control volume" annotation (Placement(transformation(extent={{-110,-88},{-90,-68}}),  iconTransformation(extent={{-110,-88},{-90,-68}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor cellBuffer(
    C=100,
    T(start=296.65, fixed=true),
    der_T(start=0))             annotation (Placement(transformation(extent={{-28,-56},{-8,-36}})));
  Modelica.Blocks.Interfaces.RealInput T_environment "Prescribed boundary temperature from weather file" annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{-80,20}})));
equation
  connect(watertowater.ports[1], pipe.port_a) annotation (Line(points={{-42,0},{-30,0}},   color={0,127,255}));
  connect(controller.y,pump. y)
    annotation (Line(points={{-57,66},{22,66},{22,14.4}},
                                                        color={0,0,127}));
  connect(pump.port_b, heatexchanger.ports[1]) annotation (Line(points={{34,0},{58,0},{58,-54},{55,-54}},
                                                                                                   color={0,127,255}));
  connect(pipe.port_b, pump.port_a) annotation (Line(points={{-10,0},{10,0}},   color={0,127,255}));
  connect(sink.ports[1], heatexchanger.ports[2]) annotation (Line(points={{68,0},{57,0},{57,-54}},   color={0,127,255}));
  connect(T_op, controller.u_s) annotation (Line(points={{-100,60},{-90,60},{-90,66},{-80,66}},
                                                                              color={0,0,127}));
  connect(pump.P, P_coolingPump) annotation (Line(points={{35.2,10.8},{64,10.8},{64,40},{100,40}},     color={0,0,127}));
  connect(heatPortCooling, cellBuffer.port) annotation (Line(points={{-100,-78},{-18,-78},{-18,-56}}, color={191,0,0}));
  connect(cellBuffer.port, heatexchanger.heatPort) annotation (Line(points={{-18,-56},{-18,-64},{46,-64}}, color={191,0,0}));
  connect(fuelcellTemperature.T, controller.u_m) annotation (Line(points={{-68,-31},{-68,54}},                                                 color={0,0,127}));
  connect(fuelcellTemperature.port, cellBuffer.port) annotation (Line(points={{-68,-52},{-68,-64},{-18,-64},{-18,-56}},           color={191,0,0}));
  connect(T_environment, watertowater.T_in) annotation (Line(points={{-100,0},{-74,0},{-74,4},{-64,4}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false), graphics={
                                          Text(
          extent={{-110,-80},{-88,-96}},
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
<p>Model parameters are adapted to a 5kW fuel cell, and fine-tuned with simulation results:</p>
<p><span style=\"font-family: Arial;\">[1] X. Xue, J. Tang, A. Smirnova et al. System level lumped-parameter dyamic modeling of PEM fuel cell. Journal of Power Sources, 133(2): 188-204. doi: 10.1016/j.jpowsour.2003.12.064.</span></p>
<p><span style=\"font-family: Arial;\">[2] J.C. Amphlett, R.F. Mann, B.A. Peppley, P.R. Roberge, A. Rodrigues. A model predicting transient responses of proton exchange membrane fuel cells. (1996) Journal of Power Sources, 61 (1-2): 183-188. doi: 10.1016/S0378-7753(96)02360-9. </span></p>
<p><span style=\"font-family: Arial;\">[3] Barbir. Fuel cell technology: Reaching towards commercialization. Engineering Materials and Processes.</span></p>
<p><span style=\"font-family: Arial;\">Initially based parameters and model for a 5kW Giner electrolyzer, based on [1] Z. Abdin, E. MacA. Gray, and C.J. Webb. Modelling and simulation of a proton exchange membrane (PEM) electrolyzer cell. International Journal of Hydrogen Energy, 40(39):13243-13257, 2015. doi:10.1016/j.ijhydene.2015.07.129. </span></p>
<p><br>Tested in the check models &quot;H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.CoolingSystem.Test_CoolingModel&quot;</p>
<p>Model created by Ali&eacute;nor Hamoir in May 2024.</p>
</html>"));
end CoolingModel;
