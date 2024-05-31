within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.OLD;
model ElyStorageController "Controller to control the electrolyzer for feeding into storage system"

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
// Gas- und WÃ¤rme-Institut Essen						  //
// and                                                                            //
// XRG Simulation GmbH (Hamburg, Germany).                                        //
//________________________________________________________________________________//

  // _____________________________________________
  //
  //          Imports and Class Hierarchy
  // _____________________________________________

  import      Modelica.Units.SI;
  extends TransiEnt.Basics.Icons.Controller;

  // _____________________________________________
  //
  //        Constants and Hidden Parameters
  // _____________________________________________

  // _____________________________________________
  //
  //             Visible Parameters
  // _____________________________________________

  parameter SI.ActivePower P_el_n=1e6 "Nominal power of the electrolyzer" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.ActivePower P_el_min=0.02*P_el_n "Maximum power of the electrolyzer" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.ActivePower P_el_max=3*P_el_n "Maximum power of the electrolyzer" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.ActivePower P_el_overload=1.5*P_el_n "Power at which the overload region of the electrolyzer begins" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.ActivePower P_el_cooldown=P_el_n "Power below which cooldown starts" annotation(Dialog(group="Fundamental Definitions"));
  parameter Modelica.Units.SI.Efficiency eta_n(
    min=0,
    max=1) = 0.75 "Nominal efficiency coefficient (min = 0, max = 1)" annotation (Dialog(group="Fundamental Definitions"));

  parameter Modelica.Units.SI.Efficiency eta_scale(
    min=0,
    max=1) = 0 "Sets a with increasing input power linear degrading efficiency coefficient (min = 0, max = 1)" annotation (Dialog(group="Fundamental Definitions"));
  parameter SI.Time t_overload=5*3600 "Maximum time in seconds that the electrolyzer can work in overload" annotation(Dialog(group="Fundamental Definitions"));
  parameter Real coolingToHeatingRatio=2 "Ratio of how much faster the electrolyzer cools down than it heats up" annotation(Dialog(group="Fundamental Definitions"));
  parameter Integer startState=1 "Initial state of the electrolyzer (1: ready to overheat, 2: working in overload, 3: cooling down)" annotation(Dialog(group="Fundamental Definitions"));

  parameter Modelica.Blocks.Types.SimpleController controllerType=Modelica.Blocks.Types.SimpleController.P "Type of controller for feed-in control" annotation (Dialog(tab="General", group="Controller"));
  parameter Real k=1 "Gain for feed-in control" annotation (Dialog(tab="General", group="Controller"));
  parameter Real Ti=0.1 "Integrator time constant for feed-in control" annotation (Dialog(tab="General", group="Controller"));
  parameter Real Td=0.1 "Derivative time constant for feed-in control" annotation (Dialog(tab="General", group="Controller"));

  parameter SI.Pressure p_minLow=20e5 "Lower limit of the target pressure in storage" annotation (Dialog(tab="General",  group="Controller"));
  parameter SI.Pressure p_maxLow=29e5 "Lower limit of the maximum pressure in storage" annotation (Dialog(tab="General", group="Controller"));
  parameter SI.Pressure p_maxHigh=30e5 "Upper limit of the maximum pressure in storage" annotation (Dialog(tab="General", group="Controller"));
  parameter SI.Pressure p_minHigh=60e5 "if valve closed and p>p_high, open valve" annotation (Dialog(tab="General", group="Controller"));
  parameter SI.Pressure p_minLow_constantDemand=50e5 "storage can be emptied via 'm_flow_hydrogenDemand_constant' up to 'p_minLow_constantDemand'" annotation (Dialog(tab="General", group="Controller"));
  parameter SI.MassFlowRate m_flow_hydrogenDemand_constant=0 "constant hydrogen demand if hydrogen is available" annotation (Dialog(tab="General", group="Controller"));

  parameter Boolean StoreAllHydrogen=false "All Hydrogen is stored before beeing fed in" annotation (Dialog(tab="General", group="Controller"));

  // _____________________________________________
  //
  //             Variable Declarations
  // _____________________________________________
  // _____________________________________________
  //
  //                 Outer Models
  // _____________________________________________

  outer TransiEnt.SimCenter simCenter;

  // _____________________________________________
  //
  //                  Interfaces
  // _____________________________________________

  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn P_el_set "Set power"    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=270,
        origin={-60,100}),iconTransformation(extent={{-20,-20},{20,20}},
        rotation=270,
        origin={-60,100})));
  TransiEnt.Basics.Interfaces.General.PressureIn p_storage "Pressure in storage" annotation (Placement(transformation(
        extent={{20,-20},{-20,20}},
        rotation=0,
        origin={100,10}),iconTransformation(extent={{120,-10},{80,30}})));
  TransiEnt.Basics.Interfaces.General.MassFlowRateIn m_flow_ely "Hydrogen mass flow out of the electrolyser" annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=180,
        origin={100,-50})));
  TransiEnt.Basics.Interfaces.General.MassFlowRateIn m_flow_feedIn "Maximum mass flow that can be fed into the natural gas system" annotation (Placement(transformation(extent={{120,50},{80,90}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_el_ely "Controlled power of the electrolyser" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={-60,-110})));

  // _____________________________________________
  //
  //           Instances of other Classes
  // _____________________________________________

  TransiEnt.Producer.Gas.Electrolyzer.Controller.FeedInStorageController feedInStorageController(
    k=k,
    p_maxLow=p_maxLow,
    p_maxHigh=p_maxHigh,
    P_el_n=P_el_n,
    P_el_max=P_el_max,
    P_el_min=P_el_min,
    eta_n=eta_n,
    eta_scale=eta_scale,
    controllerType=controllerType,
    Ti=Ti,
    Td=Td,
    startState=startState,
    redeclare model Charline = Charline) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-60,10})));

  TransiEnt.Producer.Gas.Electrolyzer.Controller.OverloadController overloadController(
    P_el_n=P_el_n,
    P_el_overload=P_el_overload,
    t_overload=t_overload,
    coolingToHeatingRatio=coolingToHeatingRatio,
    P_el_cooldown=P_el_cooldown,
    state(start=startState))
                           annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-60,-52})));

public
  replaceable model Charline = TransiEnt.Producer.Gas.Electrolyzer.Base.ElectrolyzerEfficiencyCharlineSilyzer200        constrainedby TransiEnt.Producer.Gas.Electrolyzer.Base.PartialElectrolyzerEfficiencyCharline        "Calculate the efficiency" annotation (__Dymola_choicesAllMatching=true);

equation
  // _____________________________________________
  //
  //           Characteristic Equations
  // _____________________________________________
  // _____________________________________________
  //
  //               Connect Statements
  // _____________________________________________

  connect(feedInStorageController.p, p_storage) annotation (Line(points={{-49,10},{100,10}},                        color={0,0,127}));
  connect(feedInStorageController.m_flow_feedIn, m_flow_feedIn) annotation (Line(points={{-49,18},{74,18},{74,70},{100,70}},              color={0,0,127}));
  connect(P_el_set, feedInStorageController.P_el_set) annotation (Line(points={{-60,100},{-60,21}},          color={0,0,127}));
  connect(feedInStorageController.P_el_ely, overloadController.P_el_set) annotation (Line(points={{-60,-1},{-60,-41}},         color={0,0,127}));
  connect(m_flow_ely, feedInStorageController.m_flow_bypass) annotation (Line(points={{100,-50},{74,-50},{74,2},{-49,2}},      color={0,0,127}));
  connect(P_el_ely, overloadController.P_el_ely) annotation (Line(
      points={{-60,-110},{-60,-62.8}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
  Documentation(info="<html>
<h4><span style=\"color:#008000\">1. Purpose of model</span></h4>
<p>This is a controller to control for a system with storage the electric power of the electrolyzer, the three way valve and outlet valve of the storage. it combines the FeedInStorageTWVController, FeedInStorageController and OverloadController. </p>
<h4><span style=\"color:#008000\">2. Level of detail, physical effects considered, and physical insight</span></h4>
<p>see sub models </p>
<h4><span style=\"color:#008000\">3. Limits of validity </span></h4>
<p>see sub models </p>
<h4><span style=\"color:#008000\">4. Interfaces</span></h4>
<p>P_el_set: input for the set value for the electric power </p>
<p>m_flow_feedIn: input for the possible feed-in mass flow into the natural grid etc. </p>
<p>m_flow_ely: input for the mass flow coming from the electrolyzer </p>
<p>p_storage: input for the storage pressure </p>
<p>P_el_ely: output for the limited electric power for the electrolyzer </p>
<p>splitRatio: output for the split ratio of the three way valve </p>
<p>m_flowDes_valve: output for the mass flow through the valve after the storage </p>
<h4><span style=\"color:#008000\">5. Nomenclature</span></h4>
<p>(no elements)</p>
<h4><span style=\"color:#008000\">6. Governing Equations</span></h4>
<p>The desired mass flow through valve after the storage is simply the difference between feed-in mass flow and bypass mass flow. </p>
<h4><span style=\"color:#008000\">7. Remarks for Usage</span></h4>
<p>(no remarks) </p>
<h4><span style=\"color:#008000\">8. Validation</span></h4>
<p>(no remarks) </p>
<h4><span style=\"color:#008000\">9. References</span></h4>
<p>(no remarks) </p>
<h4><span style=\"color:#008000\">10. Version History</span></h4>
<p>Model created by Carsten Bode (c.bode@tuhh.de) in April 2016<br> </p>
</html>"));
end ElyStorageController;
