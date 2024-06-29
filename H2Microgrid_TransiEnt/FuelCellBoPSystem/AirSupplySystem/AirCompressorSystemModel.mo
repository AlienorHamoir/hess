within H2Microgrid_TransiEnt.FuelCellBoPSystem.AirSupplySystem;
model AirCompressorSystemModel "Air must be fed at high pressure in FC"
parameter TransiEnt.Basics.Media.Gases.Gas_MoistAir air=TransiEnt.Basics.Media.Gases.Gas_MoistAir() "Medium model of ideal air" annotation (choicesAllMatching);
  parameter TransiEnt.Basics.Media.Gases.VLE_VDIWA_NG7_SG_O2_var medium=TransiEnt.Basics.Media.Gases.VLE_VDIWA_NG7_SG_O2_var() "Medium model of real air" annotation (choicesAllMatching);

  TransiEnt.Components.Gas.Compressor.CompressorRealGasIsothermal_L1_simple airCompressor(
    medium=medium,
    presetVariableType="m_flow",
    useMechPowerPort=true,
    m_flowInput=true,
    integrateElPower=true) annotation (Placement(transformation(extent={{-38,4},{-18,24}})));
  TransiEnt.Components.Electrical.Machines.MotorComplex
                                   motorComplex(cosphi=1, eta=0.95)
                                                          annotation (Placement(transformation(extent={{-18,-38},{2,-18}})));
  TransiEnt.Components.Boundaries.Electrical.ComplexPower.SlackBoundary
                                                   slackBoundary annotation (Placement(transformation(extent={{54,-38},{74,-18}})));
  TransiEnt.Components.Sensors.ElectricPowerComplex electricPowerComplex annotation (Placement(transformation(extent={{16,-38},{36,-18}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_pTxi boundary_pTxi(
  medium=medium,
    p_const=100000,
    T_const=296.65,
    xi_const={0,0,0,0.77,0,0.001,0,0.23,0})                                                                             annotation (Placement(transformation(extent={{-88,4},{-68,24}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_pTxi boundary_pTxi1(
    medium=medium,
    p_const=1000000,
    T_const=296.65,
    xi_const={0,0,0,0.77,0,0.001,0,0.23,0})                                                                             annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=180,
        origin={28,14})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{50,-90},{70,-70}})));
  TransiEnt.Basics.Interfaces.General.MassFlowRateIn AirMassFlowRateSetpoint annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=0,
        origin={-108,64})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_airCompressor "Active Power for the air compressor" annotation (Placement(transformation(extent={{96,-18},{116,2}})));
equation
  connect(motorComplex.epp, electricPowerComplex.epp_IN) annotation (Line(
      points={{2.1,-28.1},{2.1,-28},{16.8,-28}},
      color={28,108,200},
      thickness=0.5));
  connect(electricPowerComplex.epp_OUT, slackBoundary.epp) annotation (Line(
      points={{35.4,-28},{54,-28}},
      color={28,108,200},
      thickness=0.5));
  connect(airCompressor.mpp, motorComplex.mpp) annotation (Line(points={{-28,4},{-28,-28},{-18,-28}}, color={95,95,95}));
  connect(boundary_pTxi.gasPort, airCompressor.gasPortIn) annotation (Line(
      points={{-68,14},{-38,14}},
      color={255,255,0},
      thickness=1.5));
  connect(airCompressor.gasPortOut, boundary_pTxi1.gasPort) annotation (Line(
      points={{-18,14},{18,14}},
      color={255,255,0},
      thickness=1.5));
  connect(AirMassFlowRateSetpoint, airCompressor.m_flow_in) annotation (Line(points={{-108,64},{-36,64},{-36,25}}, color={0,0,127}));
  connect(electricPowerComplex.P, P_airCompressor) annotation (Line(
      points={{21,-19.4},{22,-19.4},{22,-8},{106,-8}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end AirCompressorSystemModel;
