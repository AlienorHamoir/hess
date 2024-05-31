within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.CoolingSystem.FluidPortCooling;
model Cooling_BuildingPump

    parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_coolant=simCenter.fluid1;
    parameter Real n_cell=20;

  TransiEnt.Basics.Interfaces.Thermal.FluidPortIn fluidPortIn(Medium=medium_coolant) annotation (Placement(transformation(extent={{-110,10},{-90,30}}), iconTransformation(extent={{90,-50},{110,-30}})));
  TransiEnt.Basics.Interfaces.Thermal.FluidPortOut fluidPortOut(Medium=medium_coolant) annotation (Placement(transformation(extent={{-110,-30},{-90,-10}}), iconTransformation(extent={{90,-90},{110,-70}})));
  TransiEnt.Components.Sensors.RealGas.MassFlowSensor massFlowSensor(medium=medium_coolant) annotation (Placement(transformation(extent={{-62,20},{-42,40}})));
  Modelica.Blocks.Interfaces.RealOutput P_el_pump "Electrical power consumed" annotation (Placement(transformation(extent={{96,64},{116,84}})));
  TransiEnt.Components.Boundaries.FluidFlow.BoundaryVLE_Txim_flow boundaryVLE_Txim_flow1(
    variable_m_flow=true,                                                                                      medium=medium_coolant,
    variable_T=false)                                                                                                                 annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={60,-20})));

  Buildings.Fluid.Sources.Boundary_pT bou1(nPorts=1, redeclare package Medium = Buildings.Media.Water "Water")
                                                                                   annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={40,66})));
  inner SimCenter simCenter annotation (Placement(transformation(extent={{-82,-92},{-62,-72}})));
  Buildings.Experimental.DHC.EnergyTransferStations.BaseClasses.Pump_m_flow pump_m_flow(
    redeclare package Medium = Buildings.Media.Water,
    m_flow_nominal=8.33*n_cell,
    dp_nominal=17.2*10e5)                                          annotation (Placement(transformation(extent={{-20,56},{0,76}})));
  TransiEnt.Components.Boundaries.FluidFlow.BoundaryVLE_pTxi boundaryVLE_pTxi(medium=medium_coolant) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={60,20})));
  Buildings.Fluid.Sources.MassFlowSource_T boundary(
    redeclare package Medium = Buildings.Media.Water,
    use_m_flow_in=true,
    nPorts=1) annotation (Placement(transformation(extent={{-66,56},{-46,76}})));
equation
  connect(fluidPortIn, massFlowSensor.gasPortIn) annotation (Line(
      points={{-100,20},{-62,20}},
      color={175,0,0},
      thickness=0.5));
  connect(boundaryVLE_Txim_flow1.fluidPortOut, fluidPortOut) annotation (Line(
      points={{50,-20},{-100,-20}},
      color={175,0,0},
      thickness=0.5));
  connect(massFlowSensor.m_flow, pump_m_flow.m_flow_in) annotation (Line(points={{-41,30},{6,30},{6,88},{-10,88},{-10,78}}, color={0,0,127}));
  connect(pump_m_flow.P, P_el_pump) annotation (Line(points={{1,75},{18,75},{18,80},{90,80},{90,74},{106,74}}, color={0,0,127}));
  connect(pump_m_flow.m_flow_actual, boundaryVLE_Txim_flow1.m_flow) annotation (Line(points={{1,71},{16,71},{16,34},{80,34},{80,-26},{72,-26}}, color={0,0,127}));
  connect(pump_m_flow.port_b, bou1.ports[1]) annotation (Line(points={{0,66},{30,66}}, color={0,127,255}));
  connect(massFlowSensor.gasPortOut, boundaryVLE_pTxi.fluidPortIn) annotation (Line(
      points={{-42,20},{50,20}},
      color={255,255,0},
      thickness=1.5));
  connect(boundary.ports[1], pump_m_flow.port_a) annotation (Line(points={{-46,66},{-20,66}}, color={0,127,255}));
  connect(boundary.m_flow_in, pump_m_flow.m_flow_actual) annotation (Line(points={{-68,74},{-76,74},{-76,50},{16,50},{16,71},{1,71}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end Cooling_BuildingPump;
