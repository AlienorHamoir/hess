within H2Microgrid_TransiEnt.FuelCellBoPSystem.AirSupplySystem;
model AirCompressorSystemModel "Air must be fed at high pressure in FC"
parameter TransiEnt.Basics.Media.Gases.Gas_MoistAir Air=TransiEnt.Basics.Media.Gases.Gas_MoistAir() "Medium model of ideal air" annotation (choicesAllMatching);
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
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_Txim_flow AirMassSource(
    medium=medium,
    variable_m_flow=false,
    m_flow_const=0.000001,
    T_const=298.15,
    xi_const={0,0,0,0.77,0,0.001,0,0.23,0}) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-78,14})));
  TransiEnt.Basics.Interfaces.General.MassFlowRateIn massFlowRateIn annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={-36,102})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_pTxi boundary_pTxi(medium=medium,
    p_const=100000,
    T_const=298.15,
    xi_const={0,0,0,0.77,0,0.001,0,0.23,0})                                             annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={32,14})));
  ClaRa.Components.TurboMachines.Compressors.CompressorGas_L1_affinity compressorGas_L1_affinity(useMechanicalPort=true) annotation (Placement(transformation(extent={{-66,-52},{-46,-32}})));
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
  connect(AirMassSource.gasPort, airCompressor.gasPortIn) annotation (Line(
      points={{-68,14},{-38,14}},
      color={255,255,0},
      thickness=1.5));
  connect(massFlowRateIn, airCompressor.m_flow_in) annotation (Line(points={{-36,102},{-36,25}}, color={0,0,127}));
  connect(airCompressor.gasPortOut, boundary_pTxi.gasPort) annotation (Line(
      points={{-18,14},{22,14}},
      color={255,255,0},
      thickness=1.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end AirCompressorSystemModel;
