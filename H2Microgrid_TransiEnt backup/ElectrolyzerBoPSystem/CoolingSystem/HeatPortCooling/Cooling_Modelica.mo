within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.CoolingSystem.HeatPortCooling;
model Cooling_Modelica

  parameter Modelica.Thermal.FluidHeatFlow.Media.Medium medium_coolant=Modelica.Thermal.FluidHeatFlow.Media.Water()
    "Cooling medium" annotation (choicesAllMatching=true);
  parameter Modelica.Units.SI.Temperature TAmb(displayUnit="degC")=293.15
    "Ambient temperature";

  Modelica.Thermal.FluidHeatFlow.Sources.Ambient
                                ambient3(
    constantAmbientTemperature=TAmb,
    medium=medium_coolant,
    constantAmbientPressure=0)
    annotation (Placement(transformation(extent={{-34,-8},{-54,12}})));
  Modelica.Thermal.FluidHeatFlow.Sources.IdealPump
                                  idealPump(
    medium=medium_coolant,
    m=0,
    T0=TAmb,
    V_flow0=0.0001272217,
    wNominal(displayUnit="rad/s") = 366.52,
    dp0(displayUnit="Pa") = 17.2*10e5)
    annotation (Placement(transformation(extent={{-24,12},{-4,-8}})));
  Modelica.Thermal.FluidHeatFlow.Components.Valve
                                 valve(
    medium=medium_coolant,
    m=0,
    T0=TAmb,
    LinearCharacteristic=false,
    y1=1,
    Kv1=1,
    kv0=0.01,
    dp0(displayUnit="Pa") = 1,
    rho0=10,
    frictionLoss=0)
    annotation (Placement(transformation(extent={{6,-8},{26,12}})));
  Modelica.Thermal.FluidHeatFlow.Components.Pipe
                                pipe1(
    medium=medium_coolant,
    T0=TAmb,
    m=0.1,
    V_flowLaminar=0.1,
    dpLaminar(displayUnit="Pa") = 0.1,
    V_flowNominal=1,
    dpNominal(displayUnit="Pa") = 1,
    h_g=0,
    T0fixed=true,
    useHeatPort=true)
    annotation (Placement(transformation(extent={{36,-8},{56,12}})));
  Modelica.Thermal.FluidHeatFlow.Sources.Ambient
                                ambient4(
    constantAmbientTemperature=TAmb,
    medium=medium_coolant,
    constantAmbientPressure=0)
    annotation (Placement(transformation(extent={{76,-8},{96,12}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor1(C=0.1, T(start=TAmb, fixed=true))
    annotation (Placement(transformation(
        origin={56,-62},
        extent={{-10,10},{10,-10}},
        rotation=90)));
  Modelica.Mechanics.Rotational.Sources.Speed speed(exact=true, useSupport=false)
    annotation (Placement(transformation(
        origin={-14,32},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  Modelica.Blocks.Sources.Ramp speedRamp(
    height=100,
    offset=0.5,
    duration=0.1,
    startTime=0.4) annotation (Placement(transformation(extent={{-52,40},{-32,60}})));
  Modelica.Blocks.Sources.Ramp valveRamp(
    height=0.5,
    offset=0.5,
    duration=0.1,
    startTime=0.9) annotation (Placement(transformation(extent={{48,42},{28,62}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=1777)   annotation (
    Placement(visible = true, transformation(origin={46,-32},   extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort annotation (Placement(transformation(extent={{36,-110},{56,-90}})));
equation
  connect(pipe1.flowPort_b,ambient4. flowPort) annotation (Line(points={{56,2},{76,2}},   color={255,0,0}));
  connect(ambient3.flowPort,idealPump. flowPort_a)
    annotation (Line(points={{-34,2},{-24,2}}, color={255,0,0}));
  connect(idealPump.flowPort_b,valve. flowPort_a)
    annotation (Line(points={{-4,2},{6,2}},    color={255,0,0}));
  connect(valve.flowPort_b,pipe1. flowPort_a) annotation (Line(points={{26,2},{36,2}},   color={255,0,0}));
  connect(speedRamp.y,speed. w_ref)
    annotation (Line(points={{-31,50},{-14,50},{-14,44}}, color={0,0,127}));
  connect(valveRamp.y,valve. y)
    annotation (Line(points={{27,52},{16,52},{16,12}},  color={0,0,127}));
  connect(speed.flange,idealPump. flange_a)
    annotation (Line(points={{-14,22},{-14,12}}));
  connect(thermalConductor.port_b, pipe1.heatPort) annotation (Line(points={{46,-22},{46,-8}}, color={191,0,0}));
  connect(heatCapacitor1.port, thermalConductor.port_a) annotation (Line(points={{46,-62},{46,-42}}, color={191,0,0}));
  connect(heatCapacitor1.port, heatPort) annotation (Line(points={{46,-62},{46,-100}}, color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end Cooling_Modelica;
