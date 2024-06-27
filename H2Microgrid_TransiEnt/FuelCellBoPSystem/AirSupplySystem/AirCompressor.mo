within H2Microgrid_TransiEnt.FuelCellBoPSystem.AirSupplySystem;
model AirCompressor "use of Modelica library"


  // Medium declaration
  replaceable package Medium = Buildings.Media.Air;
  // Parameter definition
   //Pump speed PID Controller Parameters
  parameter Modelica.Units.SI.Time tau_i=0.1 "1/tau_i for air pump system PID integrator gain ";
  parameter Real k_p=100 "gain, cooling system PID proportional control";
//   parameter Modelica.Units.SI.Time tau_d=1e-1 "tau_d, for air pump system PID derivator gain";
  parameter Real N_i=0.5 "gain of anti-windup compensation ";
//   parameter Real N_d=1 "gain, ideal derivative block ";

  Buildings.Controls.Continuous.LimPID controller(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=k_p,
    Ti=tau_i,
    Td=100,
    yMax=1,
    yMin=0,
    Ni=N_i,
    Nd=1,
    reverseActing=false)
    annotation (Placement(transformation(extent={{-66,70},{-46,90}})));
  Modelica.Blocks.Interfaces.RealOutput P_airCompressor "Electrical power consumed" annotation (Placement(transformation(extent={{96,24},{116,44}})));
  Modelica.Fluid.Sensors.MassFlowRate sen_Air_mflow(redeclare package Medium = Medium) annotation (
    Placement(visible = true, transformation(origin={40,0},      extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  TransiEnt.Basics.Interfaces.General.MassFlowRateIn AirMassFlowRateSetpoint annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=0,
        origin={-110,80})));
  Buildings.Fluid.Movers.SpeedControlled_y pump(redeclare package Medium = Medium,
    redeclare Buildings.Fluid.Movers.Data.Pumps.Wilo.Stratos25slash1to4 per,
    addPowerToMedium=false)
    annotation (Placement(transformation(extent={{-26,-12},{-2,12}})));
  Buildings.Fluid.MixingVolumes.MixingVolume FC(
    redeclare package Medium = Medium,
    m_flow_nominal=0.0001,
    V=0.1,
    nPorts=2) annotation (Placement(transformation(extent={{54,-28},{74,-48}})));
  Modelica.Blocks.Math.Gain gain(k=1)  annotation (Placement(transformation(
        extent={{-6,-6},{6,6}},
        rotation=90,
        origin={-56,34})));
  Modelica.Blocks.Math.Gain gain1(k=1) annotation (Placement(transformation(
        extent={{-4,-4},{4,4}},
        rotation=0,
        origin={-80,80})));
  Buildings.Fluid.Sources.Boundary_pT AirSource(
    redeclare package Medium = Medium,
    p=100000,
    use_T_in=false,
    T=298.15,
    nPorts=1) annotation (Placement(transformation(extent={{-92,-10},{-72,10}})));
  Buildings.Fluid.Sources.Boundary_pT AirSink(redeclare package Medium = Medium,
    p=99900,                                                                     nPorts=1) annotation (Placement(transformation(extent={{96,-10},{76,10}})));
  inner Modelica.Fluid.System system annotation (Placement(transformation(extent={{58,-90},{78,-70}})));
  Modelica.Blocks.Interfaces.RealOutput AirMassFlowToFC "Air mass flow rate from compressor to FC" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={40,-104})));
equation
  connect(pump.port_b, sen_Air_mflow.port_a) annotation (Line(points={{-2,0},{30,0}}, color={0,127,255}));
  connect(controller.y, pump.y) annotation (Line(points={{-45,80},{-14,80},{-14,14.4}}, color={0,0,127}));
  connect(pump.P, P_airCompressor) annotation (Line(points={{-0.8,10.8},{26,10.8},{26,34},{106,34}}, color={0,0,127}));
  connect(sen_Air_mflow.port_b, FC.ports[1]) annotation (Line(points={{50,0},{63,0},{63,-28}}, color={0,127,255}));
  connect(controller.u_m, gain.y) annotation (Line(points={{-56,68},{-56,40.6}}, color={0,0,127}));
  connect(gain.u, sen_Air_mflow.m_flow) annotation (Line(points={{-56,26.8},{-56,-22},{40,-22},{40,-11}}, color={0,0,127}));
  connect(AirMassFlowRateSetpoint, gain1.u) annotation (Line(points={{-110,80},{-84.8,80}}, color={0,0,127}));
  connect(gain1.y, controller.u_s) annotation (Line(points={{-75.6,80},{-68,80}}, color={0,0,127}));
  connect(AirSink.ports[1], FC.ports[2]) annotation (Line(points={{76,0},{65,0},{65,-28}}, color={0,127,255}));
  connect(AirSource.ports[1], pump.port_a) annotation (Line(points={{-72,0},{-26,0}}, color={0,127,255}));
  connect(sen_Air_mflow.m_flow, AirMassFlowToFC) annotation (Line(points={{40,-11},{40,-104}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end AirCompressor;
