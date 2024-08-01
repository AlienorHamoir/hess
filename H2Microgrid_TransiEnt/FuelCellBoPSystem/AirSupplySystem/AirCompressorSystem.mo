within H2Microgrid_TransiEnt.FuelCellBoPSystem.AirSupplySystem;
model AirCompressorSystem "Air compressor system model"
  parameter TransiEnt.Basics.Media.Gases.Gas_MoistAir air=TransiEnt.Basics.Media.Gases.Gas_MoistAir() "Medium model of ideal air" annotation (choicesAllMatching);
  parameter TransiEnt.Basics.Media.Gases.VLE_VDIWA_NG7_SG_O2_var medium=TransiEnt.Basics.Media.Gases.VLE_VDIWA_NG7_SG_O2_var() "Medium model of real air" annotation (choicesAllMatching); // in this model, we use real air (on the contrary in fuel cell model, ideal air is used)
  parameter Modelica.Units.SI.Efficiency eta_mech_compressor(
    min=0,
    max=1)=0.95 "Compressor mechanical efficiency coefficient (min = 0, max = 1)" annotation (Dialog(tab="General", group="Compressor"));
  parameter Modelica.Units.SI.Efficiency eta_el_compressor(
    min=0,
    max=1)=0.95 "Compressor motor electrical efficiency coefficient (min = 0, max = 1)" annotation (Dialog(tab="General", group="Compressor"));
  TransiEnt.Components.Gas.Compressor.CompressorRealGasIsothermal_L1_simple airCompressor(
    medium=medium,
    eta_mech=eta_mech_compressor,
    eta_el=eta_el_compressor,
    presetVariableType="m_flow",
    useMechPowerPort=true,
    m_flowInput=true,
    integrateElPower=true,
    allow_reverseFlow=false)
                           annotation (Placement(transformation(extent={{-38,4},{-18,24}})));
  TransiEnt.Components.Electrical.Machines.MotorComplex
                                   motorComplex(cosphi=1, eta=0.95)
                                                          annotation (Placement(transformation(extent={{-18,-50},{2,-30}})));
  TransiEnt.Components.Boundaries.Electrical.ComplexPower.SlackBoundary
                                                   slackBoundary annotation (Placement(transformation(extent={{54,-50},{74,-30}})));
  TransiEnt.Components.Sensors.ElectricPowerComplex electricPowerComplex annotation (Placement(transformation(extent={{16,-50},{36,-30}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_pTxi boundary_pTxi(
  medium=medium,
    p_const=100000,
    T_const=296.65,
    xi_const={0,0,0,0.77,0,0.02,0,0.21,0})                                                                              annotation (Placement(transformation(extent={{-88,4},{-68,24}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_pTxi boundary_pTxi1(
    medium=medium,
    p_const=300000,
    T_const=296.65,
    xi_const={0,0,0,0.77,0,0.02,0,0.21,0})                                                                              annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=180,
        origin={28,14})));
  TransiEnt.Basics.Interfaces.General.MassFlowRateIn AirMassFlowRateSetpoint annotation (Placement(transformation(
        extent={{-16,-16},{16,16}},
        rotation=0,
        origin={-100,60}), iconTransformation(extent={{-16,-16},{16,16}}, origin={-100,60})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_airCompressor "Active Power for the air compressor" annotation (Placement(transformation(extent={{90,-30},{110,-10}})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{-86,-94},{-66,-74}})));
equation
  connect(motorComplex.epp, electricPowerComplex.epp_IN) annotation (Line(
      points={{2.1,-40.1},{2.1,-40},{16.8,-40}},
      color={28,108,200},
      thickness=0.5));
  connect(electricPowerComplex.epp_OUT, slackBoundary.epp) annotation (Line(
      points={{35.4,-40},{54,-40}},
      color={28,108,200},
      thickness=0.5));
  connect(airCompressor.mpp, motorComplex.mpp) annotation (Line(points={{-28,4},{-28,-40},{-18,-40}}, color={95,95,95}));
  connect(boundary_pTxi.gasPort, airCompressor.gasPortIn) annotation (Line(
      points={{-68,14},{-38,14}},
      color={255,255,0},
      thickness=1.5));
  connect(airCompressor.gasPortOut, boundary_pTxi1.gasPort) annotation (Line(
      points={{-18,14},{18,14}},
      color={255,255,0},
      thickness=1.5));
  connect(electricPowerComplex.P, P_airCompressor) annotation (Line(
      points={{21,-31.4},{21,-20},{100,-20}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(AirMassFlowRateSetpoint, airCompressor.m_flow_in) annotation (Line(points={{-100,60},{-36,60},{-36,25}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<h4>1. Purpose of model</h4>
<p>Air compressor model for FC applications. </p>
<p>Air must be fed at high pressure and at an adequate mass flow, computed by the OER controller, for the electrochemical reaction to take place. </p>
<p>It was decided air mus be fed at around 10 bar in FC, which is equal to the input pressure of hydrogen for EH Group 20-40 kW EH-Trace PEM fuel cells</p>
<p>System consists of a real gas isothermal simple compressor model, and its flange motor. Components were chosen in TransiEnt library, for compatibility with the rest of the model.</p>
<p>A real gas was used to modelize air, with adapted composition: Ni = 0.77, O2 = 0.22, H2O = 0.001.</p>
<p>2. Level of detail, physical effects considered, and physical insight</p>
<p>(Purely technical component without physical modeling.)</p>
<h4>3. Limits of validity </h4>
<p>(Purely technical component without physical modeling.)</p>
<h4>4. Interfaces</h4>
<p>(no remarks)</p>
<h4>5. Nomenclature</h4>
<p>(no elements)</p>
<h4>6. Governing Equations</h4>
<p>(no equations)</p>
<h4>7. Remarks for Usage</h4>
<p>(no remarks)</p>
<h4>8. Validation</h4>
<p>Tested in the check models &quot;H2Microgrid_TransiEnt.FuelCellBoPSystem.AirSupplySystem.TestAirCompressor&quot;</p>
<h4>9. References</h4>
<p>(no remarks)</p>
<h4>10. Version History</h4>
<p>Model created by Ali&eacute;nor Hamoir in June 2024.</p>
</html>"));
end AirCompressorSystem;
