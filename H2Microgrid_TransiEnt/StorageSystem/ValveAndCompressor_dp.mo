within H2Microgrid_TransiEnt.StorageSystem;
model ValveAndCompressor_dp

//________________________________________________________________________________//
// Component of the TransiEnt Library, version: 2.0.3                             //
//                                                                                //
// Licensed by Hamburg University of Technology under the 3-BSD-clause.           //
// Copyright 2021, Hamburg University of Technology.                              //
//________________________________________________________________________________//
//                                                                                //
// TransiEnt.EE, ResiliEntEE, IntegraNet and IntegraNet II are research projects  //
// supported by the German Federal Ministry of Economics and Energy               //
// (FKZ 03ET4003, 03ET4048, 0324027 and 03EI1008).                                //
// The TransiEnt Library research team consists of the following project partners://
// Institute of Engineering Thermodynamics (Hamburg University of Technology),    //
// Institute of Energy Systems (Hamburg University of Technology),                //
// Institute of Electrical Power and Energy Technology                            //
// (Hamburg University of Technology)                                             //
// Fraunhofer Institute for Environmental, Safety, and Energy Technology UMSICHT, //
// Gas- und WÃ¤rme-Institut Essen                                                  //
// and                                                                            //
// XRG Simulation GmbH (Hamburg, Germany).                                        //
//________________________________________________________________________________//

  // _____________________________________________
  //
  //          Imports and Class Hierarchy
  // _____________________________________________

  import      Modelica.Units.SI;
  extends TransiEnt.Basics.Icons.Model;

  // _____________________________________________
  //
  //        Constants and Hidden Parameters
  // _____________________________________________

  // _____________________________________________
  //
  //             Visible Parameters
  // _____________________________________________

  replaceable parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3 "Hydrogen to be used" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.Volume volumeSplit = 0.1 "Volume of the split" annotation(Dialog(group="Split and Junction"));
  parameter SI.Volume volumeJunction = 0.1 "Volume of the junction" annotation(Dialog(group="Split and Junction"));

  parameter SI.Pressure p_startSplit = 1e5 "Start pressure in the split" annotation(Dialog(group="Initial Values"));
  parameter SI.Pressure p_startJunction = 1e5 "Start pressure in the junction" annotation(Dialog(group="Initial Values"));
  parameter SI.MassFraction xi_startSplit[medium.nc-1] = medium.xi_default "Start composition in the split" annotation(Dialog(group="Initial Values"));
  parameter SI.MassFraction xi_startJunction[medium.nc-1] = medium.xi_default "Start composition in the junction" annotation(Dialog(group="Initial Values"));
  parameter SI.Temperature T_startSplit = simCenter.T_ground "Start temperature in the split" annotation(Dialog(group="Initial Values"));
  parameter SI.Temperature T_startJunction = simCenter.T_ground "Start temperature in the Junction" annotation(Dialog(group="Initial Values"));
  parameter SI.SpecificEnthalpy h_startSplit = TILMedia.Internals.VLEFluidConfigurations.FullyMixtureCompatible.VLEFluidFunctions.specificEnthalpy_pTxi(medium, p_startSplit, T_startSplit, xi_startSplit) "Start specific enthalpy in the split" annotation(Dialog(group="Initial Values"));
  parameter SI.SpecificEnthalpy h_startJunction = TILMedia.Internals.VLEFluidConfigurations.FullyMixtureCompatible.VLEFluidFunctions.specificEnthalpy_pTxi(medium, p_startJunction, T_startJunction, xi_startJunction) "Start specific enthalpy in the Junction" annotation(Dialog(group="Initial Values"));
  parameter Integer initOptionSplit=0 "Type of initialisation" annotation (Dialog(tab="Initialisation"), choices(
      choice=0 "Use guess values",
      choice=1 "Steady state",
      choice=201 "Steady pressure",
      choice=202 "Steady enthalpy",
      choice=208 "Steady pressure and enthalpy",
      choice=210 "Steady density"));
  parameter Integer initOptionJunction=0 "Type of initialisation" annotation (Dialog(tab="Initialisation"), choices(
      choice=0 "Use guess values",
      choice=1 "Steady state",
      choice=201 "Steady pressure",
      choice=202 "Steady enthalpy",
      choice=208 "Steady pressure and enthalpy",
      choice=210 "Steady density"));
  parameter Boolean useFluidModelsForSummary=false "True, if fluid models shall be used for the summary" annotation(Dialog(tab="Summary"));

  // _____________________________________________
  //
  //                 Outer Models
  // _____________________________________________

  outer TransiEnt.SimCenter simCenter;
  replaceable model Compressor = TransiEnt.Components.Gas.Compressor.CompressorRealGasIsothermal_L1_simple constrainedby TransiEnt.Components.Gas.Compressor.Base.PartialCompressorRealGas_L1_simple                   "Compressor model, change eta_is for normal compressor" annotation(Dialog(group="Compressor"), choicesAllMatching);

  // _____________________________________________
  //
  //             Variable Declarations
  // _____________________________________________

  // _____________________________________________
  //
  //                  Interfaces
  // _____________________________________________

public
  TransiEnt.Basics.Interfaces.Gas.RealGasPortOut gasPortOut(Medium=medium) annotation (Placement(transformation(rotation=0, extent={{90,-10},{110,10}})));
  TransiEnt.Basics.Interfaces.Gas.RealGasPortIn gasPortIn(Medium=medium) annotation (Placement(transformation(rotation=0, extent={{-110,-10},{-90,10}})));
  TransiEnt.Basics.Interfaces.General.PressureDifferenceIn dp_desired "Desired pressure difference" annotation (Placement(transformation(
        rotation=270,
        extent={{-10,-10},{10,10}},
        origin={0,90})));
  // _____________________________________________
  //
  //           Instances of other Classes
  // _____________________________________________

protected
  model ValveRecord
    extends TransiEnt.Basics.Icons.Record;
    parameter Boolean showSummary=true;

    input SI.MassFlowRate m_flow if showSummary "Mass flow rate";
    input SI.PressureDifference Delta_p if showSummary "Pressure difference";
    input SI.Temperature T_in if showSummary "Inlet temperature";
    input SI.Temperature T_out if showSummary "Outlet temperature";
    input SI.Pressure p_in if showSummary "Inlet pressure";
    input SI.Pressure p_out if showSummary "Outlet pressure";
    input SI.SpecificEnthalpy h_in if showSummary "Inlet specific enthalpy";
    input SI.SpecificEnthalpy h_out if showSummary "Outlet specific enthalpy";
  end ValveRecord;

  model CompressorRecord
    extends TransiEnt.Basics.Icons.Record;
    input SI.MassFlowRate m_flow "Mass flow rate";
    input SI.VolumeFlowRate V_flow "Volume flow rate";
    input SI.Power P_hyd "Hydraulic power";
    input SI.Power P_shaft "Shaft power";
    input SI.Power P_el "Electric power";
    input SI.Work W_el "Electric work";
    input Real Pi "Pressure ratio";
    input SI.PressureDifference Delta_p "Pressure difference";
    input Real eta "Isentropic efficiency";
    input SI.Temperature T_in "Inlet temperature";
    input SI.Temperature T_out "Outlet temperature";
    input SI.Pressure p_in "Inlet pressure";
    input SI.Pressure p_out "Outlet pressure";
    input SI.SpecificEnthalpy h_in "Inlet specific enthalpy";
    input SI.SpecificEnthalpy h_out "Outlet specific enthalpy";
  end CompressorRecord;

  inner model Summary
    extends TransiEnt.Basics.Icons.Record;
    TransiEnt.Basics.Records.FlangeRealGas gasPortIn;
    TransiEnt.Basics.Records.FlangeRealGas gasPortOut;
    CompressorRecord compressor;
    ValveRecord valve;
    TransiEnt.Basics.Records.Costs costs;
  end Summary;

//   inner Summary summary(
//     gasPortIn(
//       mediumModel=medium,
//       xi=split.summary.gasPort1.xi,
//       x=split.summary.gasPort1.x,
//       m_flow=split.summary.gasPort1.m_flow,
//       T=split.summary.gasPort1.T,
//       p=split.summary.gasPort1.p,
//       h=split.summary.gasPort1.h,
//       rho=split.summary.gasPort1.rho),
//     gasPortOut(
//       mediumModel=medium,
//       xi=junction.summary.gasPort3.xi,
//       x=junction.summary.gasPort3.x,
//       m_flow=junction.summary.gasPort3.m_flow,
//       T=junction.summary.gasPort3.T,
//       p=junction.summary.gasPort3.p,
//       h=junction.summary.gasPort3.h,
//       rho=junction.summary.gasPort3.rho),
//     compressor(
//       m_flow=compressor.summary.gasPortIn.m_flow,
//       V_flow=compressor.summary.outline.V_flow,
//       P_hyd=compressor.summary.outline.P_hyd,
//       P_shaft=compressor.summary.outline.P_shaft,
//       P_el=compressor.summary.outline.P_el,
//       W_el=compressor.summary.outline.W_el,
//       Pi=compressor.summary.outline.Pi,
//       Delta_p=compressor.summary.outline.Delta_p,
//       eta=compressor.summary.outline.eta,
//       T_in=compressor.summary.gasPortIn.T,
//       T_out=compressor.summary.gasPortOut.T,
//       p_in=compressor.summary.gasPortIn.p,
//       p_out=compressor.summary.gasPortOut.p,
//       h_in=compressor.summary.gasPortIn.h,
//       h_out=compressor.summary.gasPortOut.h),
//     valve(
//       m_flow=valve.summary.gasPortIn.m_flow,
//       showSummary=useFluidModelsForSummary,
//       Delta_p=valve.gasPortIn.p - gasPortOut.p,
//       T_in=valve.summary.gasPortIn.T,
//       T_out=valve.summary.gasPortOut.T,
//       p_in=valve.summary.gasPortIn.p,
//       p_out=valve.summary.gasPortOut.p,
//       h_in=valve.summary.gasPortIn.h,
//       h_out=valve.summary.gasPortOut.h),
//     costs(
//       costs=compressor.summary.costs.costs,
//       investCosts=compressor.summary.costs.investCosts,
//       demandCosts=compressor.summary.costs.demandCosts,
//       oMCosts=compressor.summary.costs.oMCosts,
//       otherCosts=compressor.summary.costs.otherCosts,
//       revenues=compressor.summary.costs.revenues)) annotation (Placement(transformation(extent={{-100,-100},{-80,-80}})));

public
  ControllerValveAndCompressor_dp controllerValveAndCompressor_dp annotation (Placement(transformation(extent={{-10,42},{10,62}})));
  TransiEnt.Components.Gas.Compressor.CompressorRealGasIsentropicEff_L1_simple compressorRealGasIsentropicEff_L1_simple(
    medium=medium,
    eta_mech=0.9,
    eta_el=0.9,
    presetVariableType="dp",
    useMechPowerPort=true,
    use_Delta_p_input=true,
    P_el_n=60,
    eta_is=0.7) annotation (Placement(transformation(extent={{-8,-10},{12,10}})));
  TransiEnt.Basics.Interfaces.General.MechanicalPowerPort mpp2 annotation (Placement(transformation(extent={{58,-102},{78,-82}})));
protected
  TransiEnt.Components.Gas.VolumesValvesFittings.Fittings.RealGasJunction_L2 junction(
    medium=medium,
    volume=volumeJunction,
    initOption=initOptionJunction,
    p_start=p_startJunction,
    h_start=h_startJunction,
    xi_start=xi_startJunction,
    T_start=T_startJunction) annotation (Placement(transformation(extent={{28,-10},{48,10}})));
  TransiEnt.Components.Gas.VolumesValvesFittings.Fittings.RealGasJunction_L2 split(
    medium=medium,
    volume=volumeSplit,
    initOption=0,
    p_start=p_startSplit,
    h_start=h_startSplit,
    xi_start=xi_startSplit,
    T_start=T_startSplit) annotation (Placement(transformation(extent={{-42,-10},{-22,10}})));

  TransiEnt.Components.Gas.VolumesValvesFittings.Valves.ValveDesiredMassFlow valve(medium=medium, useFluidModelsForSummary=useFluidModelsForSummary) annotation (Placement(transformation(extent={{-6,-58},{14,-46}})));
  TransiEnt.Components.Sensors.RealGas.MassFlowSensor massFlowSensoreBefore(medium=medium) annotation (Placement(transformation(extent={{-76,0},{-56,20}})));

equation
  // _____________________________________________
  //
  //           Characteristic Equations
  // _____________________________________________

  connect(junction.gasPort2, valve.gasPortOut) annotation (Line(
      points={{38,-10},{38,-52.8571},{14,-52.8571}},
      color={255,255,0},
      thickness=1.5));
  connect(valve.gasPortIn, split.gasPort2) annotation (Line(
      points={{-6,-52.8571},{-32,-52.8571},{-32,-10}},
      color={255,255,0},
      thickness=1.5));
  connect(split.gasPort1, massFlowSensoreBefore.gasPortOut) annotation (Line(
      points={{-42,0},{-56,0}},
      color={255,255,0},
      thickness=1.5));
  connect(junction.gasPort3, gasPortOut) annotation (Line(
      points={{48,0},{100,0}},
      color={255,255,0},
      thickness=1.5));
  connect(gasPortIn, massFlowSensoreBefore.gasPortIn) annotation (Line(
      points={{-100,0},{-76,0}},
      color={255,255,0},
      thickness=1.5));
  connect(dp_desired, controllerValveAndCompressor_dp.dp_desired) annotation (Line(points={{0,90},{0,62}},  color={0,0,127}));
  connect(massFlowSensoreBefore.m_flow, controllerValveAndCompressor_dp.m_flow) annotation (Line(points={{-55,10},{-44,10},{-44,52},{-10,52}}, color={0,0,127}));
  connect(controllerValveAndCompressor_dp.m_flow_valve, valve.m_flowDes) annotation (Line(points={{-4,41},{-4,14},{-14,14},{-14,-42},{-6,-42},{-6,-47.7143}},   color={0,0,127}));
  connect(split.gasPort3, compressorRealGasIsentropicEff_L1_simple.gasPortIn) annotation (Line(
      points={{-22,0},{-8,0}},
      color={255,255,0},
      thickness=1.5));
  connect(compressorRealGasIsentropicEff_L1_simple.gasPortOut, junction.gasPort1) annotation (Line(
      points={{12,0},{28,0}},
      color={255,255,0},
      thickness=1.5));
  connect(compressorRealGasIsentropicEff_L1_simple.dp_in, controllerValveAndCompressor_dp.dp_comp) annotation (Line(points={{10,11},{10,32},{4,32},{4,41}},
                                                                                                                                                          color={0,0,127}));
  connect(compressorRealGasIsentropicEff_L1_simple.mpp, mpp2) annotation (Line(points={{2,-10},{2,-16},{68,-16},{68,-92}}, color={95,95,95}));
  annotation (Diagram(coordinateSystem(extent={{-100,-100},{100,100}})), Icon(coordinateSystem(extent={{-100,-100},{100,100}}), graphics={
        Ellipse(extent={{-40,80},{40,0}}, lineColor={0,0,0}),
        Line(points={{-40,-20},{-40,-60},{40,-20},{40,-60},{-40,-20}}, color={0,0,0}),
        Line(points={{-30,66},{40,46}}, color={0,0,0}),
        Line(points={{-30,14},{40,34}}, color={0,0,0}),
        Line(
          points={{-40,40},{-70,40},{-70,0},{-100,0},{-70,0},{-70,-40},{-40,-40}},
          color={255,255,0},
          thickness=0.5),
        Line(
          points={{40,40},{70,40},{70,0},{100,0},{70,0},{70,-40},{40,-40}},
          color={255,255,0},
          thickness=0.5)}),
          Documentation(info="<html>
<h4><span style=\"color: #008000\">1. Purpose of model</span></h4>
<p>This model combines a valve with a compressor delivering a certain pressure difference and a certain mass flow rate. Depending on the sign of the pressure difference the compressor or the valve is used. The type of compressor can be chosen.</p>
<h4><span style=\"color: #008000\">2. Level of detail, physical effects considered, and physical insight</span></h4>
<p>When the pressure difference changes sign, the change of the device is executed without any time delay or time dependent behaviour. </p>
<p>From initial model, a motor port has been added for assesment of compressor power consumption. The medium has been changed to pure hydrogen.</p>
<h4><span style=\"color: #008000\">3. Limits of validity </span></h4>
<p>Only valid for negligible time delays and time dependencies in turn-on and shut-down behaviour of the components. </p>
<h4><span style=\"color: #008000\">4. Interfaces</span></h4>
<p>gasPortIn: real gas inlet </p>
<p>gasPortOut: real gas outlet </p>
<p>dp_desired: input for desired pressure difference in [Pa]</p>
<h4><span style=\"color: #008000\">5. Nomenclature</span></h4>
<p>(no elements)</p>
<h4><span style=\"color: #008000\">6. Governing Equations</span></h4>
<p>(no remarks) </p>
<h4><span style=\"color: #008000\">7. Remarks for Usage</span></h4>
<p>Chattering can occur if the pressure difference is around zero.</p>
<h4><span style=\"color: #008000\">8. Validation</span></h4>
<p>(no remarks) </p>
<h4><span style=\"color: #008000\">9. References</span></h4>
<p>Source: TransiEnt library</p>
<h4><span style=\"color: #008000\">10. Version History</span></h4>
<p>Model created by Carsten Bode (c.bode@tuhh.de) on Feb 10 2017</p>
<p>Modified by Lisa Andresen (andresen@tuhh.de) on Feb 16 2017 (changed output from m_flowDesired to dp_desired)</p>
<p>Model revised by Carsten Bode (c.bode@tuhh.de) in Apr 2018 (fixed for update to ClaRa 1.3.0)</p>
<p>Model revised by Ali&eacute;nor Hamoir in June 2024.</p>
</html>"));
end ValveAndCompressor_dp;
