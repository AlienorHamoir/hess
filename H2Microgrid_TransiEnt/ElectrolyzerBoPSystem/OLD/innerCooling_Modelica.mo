within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.OLD;
model innerCooling_Modelica
  import Modelica;

  parameter Modelica.Thermal.FluidHeatFlow.Media.Medium outerMedium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC()
    "Outer medium" annotation (choicesAllMatching=true);
  parameter Modelica.Thermal.FluidHeatFlow.Media.Medium innerMedium=Modelica.Thermal.FluidHeatFlow.Media.Water()
    "Inner medium" annotation (choicesAllMatching=true);

  parameter Modelica.Units.SI.Temperature T_amb=23 + 273.15 "K, ambient temperature";
  parameter Modelica.Units.SI.Temperature T_cool_circuit=19 + 273.15 "K, cooling circuit temperature";

  parameter Integer n_cells=Specification.n_cells "Number of electrolysis cells in series";

  replaceable H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.Specifications.GinerELX5kW Specification constrainedby H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.Specifications.Base5kWElectrolyzerL2Specification;

  output Modelica.Units.SI.TemperatureDifference dTCooler=innerPipe.T_q-outerPipe.T_q
    "Cooler's temperature increase between inner and outer pipes";
  output Modelica.Units.SI.TemperatureDifference dTouterCoolant=outerPipe.dT
    "Outer coolant's temperature increase";

  Modelica.Thermal.FluidHeatFlow.Sources.Ambient
                                ambient1(
    useTemperatureInput=false,
    constantAmbientTemperature=T_cool_circuit,
    medium=outerMedium,
    constantAmbientPressure=150000)
    annotation (Placement(transformation(extent={{-34,52},{-54,72}})));
  Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow
                                   outerPump(
    medium=outerMedium,
    m=0,
    T0=T_cool_circuit,
    useVolumeFlowInput=false,
    constantVolumeFlow=1)
    annotation (Placement(transformation(extent={{-14,52},{6,72}})));
  Modelica.Thermal.FluidHeatFlow.Sources.Ambient
                                ambient2(
    useTemperatureInput=false,
    constantAmbientTemperature=T_cool_circuit,
    medium=outerMedium,
    constantAmbientPressure=150000)
    annotation (Placement(transformation(extent={{66,52},{86,72}})));
  Modelica.Thermal.FluidHeatFlow.Sources.AbsolutePressure
                                         absolutePressure(p=450000, medium=innerMedium)
    annotation (Placement(transformation(extent={{66,-68},{86,-48}})));
  Modelica.Blocks.Sources.Constant outerVolumeFlow(k=0.0010603)
    annotation (Placement(transformation(extent={{-34,72},{-14,92}})));
  Modelica.Blocks.Sources.Constant outerGc(k=100)
    annotation (Placement(transformation(extent={{-14,22},{6,42}})));
  Modelica.Thermal.FluidHeatFlow.Components.Pipe
                                outerPipe(
    medium=outerMedium,
    m=0.1,
    T0=T_cool_circuit,
    tapT=T_cool_circuit,
    V_flowLaminar=0.1,
    dpLaminar(displayUnit="Pa") = 0,
    V_flowNominal=1,
    dpNominal(displayUnit="Pa") = 0,
    h_g=0,
    T0fixed=true,
    useHeatPort=true)
    annotation (Placement(transformation(extent={{26,52},{46,72}})));
  Modelica.Thermal.FluidHeatFlow.Components.Pipe
                                innerPipe(
    medium=innerMedium,
    m=1,
    T0=T_amb,
    V_flowLaminar=1*10e-6*n_cells,
    dpLaminar(displayUnit="Pa") = 0,
    V_flowNominal=8.33*10e-6*n_cells,
    dpNominal(displayUnit="Pa") = 0,
    h_g=0,
    T0fixed=true,
    useHeatPort=true)
    annotation (Placement(transformation(extent={{26,-36},{46,-56}})));
  Modelica.Thermal.HeatTransfer.Components.Convection outerConvection
    annotation (Placement(transformation(
        origin={36,32},
        extent={{10,10},{-10,-10}},
        rotation=270)));
  Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_a flowPort_in(medium=innerMedium)                                  annotation (Placement(transformation(extent={{-24,-110},{-4,-90}})));
  Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_b flowPort_out(medium=innerMedium)                                  annotation (Placement(transformation(extent={{56,-108},{76,-88}})));
  Modelica.Mechanics.Rotational.Sources.Speed speed(exact=true, useSupport=false)
    annotation (Placement(transformation(
        origin={-48,-68},
        extent={{-10,-10},{10,10}},
        rotation=0)));
  Modelica.Blocks.Sources.Ramp speedRamp(
    height=3500,
    offset=0.5,
    duration=0.1,
    startTime=0.1) annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-80,-68})));
  Modelica.Thermal.FluidHeatFlow.Sources.IdealPump
                                  idealPump(
    medium=innerMedium,
    m=1,
    T0=T_amb,
    T0fixed=false,
    V_flow(start=6.67*10e-6*n_cells),
    V_flow0=0.0001272217,
    wNominal(displayUnit="rad/s") = 366.52,
    dp0(displayUnit="Pa") = 17.2*10e5)
    annotation (Placement(transformation(extent={{-10,10},{10,-10}},
        rotation=90,
        origin={-14,-68})));
  TransiEnt.Producer.Heat.Heat2Heat.Indirect_HEX_const_T_out_L1 indirect_HEX_const_T_out_L1_1(
    T_start=T_amb,
    Q_max=5000,
    m_flow_max=16.67*10e-3*n_cells,
    m_flow_min=6.67*10e-3*n_cells)
                                  annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={-174,-12})));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor annotation (Placement(transformation(
        extent={{-6,-6},{6,6}},
        rotation=-90,
        origin={36,2})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor annotation (Placement(transformation(
        extent={{-6,-6},{6,6}},
        rotation=180,
        origin={12,-12})));
  inner SimCenter simCenter annotation (Placement(transformation(extent={{-84,2},{-64,22}})));
equation
  connect(ambient1.flowPort,outerPump. flowPort_a)
    annotation (Line(points={{-34,62},{-14,62}}, color={255,0,0}));
  connect(innerPipe.flowPort_b,absolutePressure. flowPort) annotation (Line(
        points={{46,-46},{60,-46},{60,-58},{66,-58}},
                                             color={255,0,0}));
  connect(outerPump.flowPort_b,outerPipe. flowPort_a)
    annotation (Line(points={{6,62},{26,62}},  color={255,0,0}));
  connect(outerPipe.flowPort_b,ambient2. flowPort)
    annotation (Line(points={{46,62},{66,62}}, color={255,0,0}));
  connect(outerPipe.heatPort,outerConvection. fluid)
    annotation (Line(points={{36,52},{36,42}},         color={191,0,0}));
  connect(outerGc.y,outerConvection. Gc)
    annotation (Line(points={{7,32},{26,32}},  color={0,0,127}));
  connect(outerVolumeFlow.y,outerPump. volumeFlow) annotation (Line(
      points={{-13,82},{-4,82},{-4,72}},   color={0,0,127}));
  connect(flowPort_out, flowPort_out) annotation (Line(points={{66,-98},{66,-98}},   color={255,0,0}));
  connect(speedRamp.y,speed. w_ref)
    annotation (Line(points={{-69,-68},{-60,-68}},        color={0,0,127}));
  connect(speed.flange,idealPump. flange_a)
    annotation (Line(points={{-38,-68},{-24,-68}}));
  connect(idealPump.flowPort_b, innerPipe.flowPort_a) annotation (Line(points={{-14,-58},{-14,-46},{26,-46}},color={255,0,0}));
  connect(absolutePressure.flowPort, flowPort_out) annotation (Line(points={{66,-58},{62,-58},{62,-84},{66,-84},{66,-98}}, color={255,0,0}));
  connect(innerPipe.heatPort, heatFlowSensor.port_b) annotation (Line(points={{36,-36},{36,-4}},  color={191,0,0}));
  connect(temperatureSensor.port, innerPipe.heatPort) annotation (Line(points={{18,-12},{36,-12},{36,-36}},
                                                                                                   color={191,0,0}));
  connect(flowPort_in, idealPump.flowPort_a) annotation (Line(points={{-14,-100},{-14,-78}}, color={255,0,0}));
  connect(heatFlowSensor.port_a, outerConvection.solid) annotation (Line(points={{36,8},{36,22}}, color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end innerCooling_Modelica;
