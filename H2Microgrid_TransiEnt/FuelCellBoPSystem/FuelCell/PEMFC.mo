within H2Microgrid_TransiEnt.FuelCellBoPSystem.FuelCell;
model PEMFC "Model of PEM Fuel Cell stack"

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
// Gas- und WÃ¤rme-Institut Essen                                                 //
// and                                                                            //
// XRG Simulation GmbH (Hamburg, Germany).                                        //
//________________________________________________________________________________//

  // _____________________________________________
  //
  //          Imports and Class Hierarchy
  // _____________________________________________
  extends TransiEnt.Basics.Icons.Model;
  import TransiEnt;


  // _____________________________________________
  //
  //                 Outer Models
  // _____________________________________________

   outer TransiEnt.SimCenter simCenter;

  // _____________________________________________
  //
  //             Visible Parameters
  // _____________________________________________

  parameter Integer no_Cells = 35 "Number of cells connected in series";

  parameter Real lambda=14 "constant humidity";

  parameter Modelica.Units.SI.Thickness t_mem = 1.275e-4 "PE membrane thickness - ref. Nafion 117";

//   parameter Real Ri = 0.15e-4 "Internal resistance according to Barbir in Ohm - here we use Ri = t_mem / mem_conductivity * A_cell or Ri(T_stack;  I) ";

  parameter Real z = 2 "Quantity of transfered electrons";

  parameter Modelica.Units.SI.Area A_cell=0.0232 "Area of one cell";

  parameter Modelica.Units.SI.Pressure p_Anode=2.4318e5 "Pressure at the anode";

  parameter Modelica.Units.SI.Pressure p_Kathode=2.4318e5 "Pressure at the cathode";

  parameter Modelica.Units.SI.Pressure p_Amb=1e5 "Pressure at the cathode";

  parameter TransiEnt.Basics.Media.Gases.Gas_VDIWA_SG7_var Syngas=TransiEnt.Basics.Media.Gases.Gas_VDIWA_SG7_var() "Medium model of Syngas" annotation (choicesAllMatching);

  parameter TransiEnt.Basics.Media.Gases.Gas_MoistAir Air=TransiEnt.Basics.Media.Gases.Gas_MoistAir() "Medium model of air - Moist air is used which consists of N2, H2O and O2. O2 is used in FC" annotation (choicesAllMatching);

  parameter Real cp_tab[3,3] = [28.91404, -0.00084, 2.01e-6; 25.84512, 0.012987, -3.9e-6; 30.62644, 0.009621, 1.18e-6] "Empiric parameter for calculating the Gibbs free energy according to Barbir";

  parameter Modelica.Units.SI.Current I_shutdown=20 "If load controller requests currents below this value, stack will shut down";

  parameter Modelica.Units.SI.Mass m=43 "Mass of the stack";

  parameter Modelica.Units.SI.SpecificHeatCapacity cp=35000/m "Specific heat capacity of the stack";

  parameter Modelica.Units.SI.Temperature T_nom = 40 + 273.15 "Temperature in nominal point";

  parameter Modelica.Units.SI.Temperature T_amb = 23.5 + 273.15 "Ambient temperature";

  parameter Modelica.Units.SI.Temperature T_stack_max = 40 + 273.15 "Maximum stack temperature";

  parameter Modelica.Units.SI.Temperature T_cool_set = 30 + 273.15 "Cooling trigger point";

  parameter Modelica.Units.SI.ThermalConductance ka=17 "Thermal conductance of the stack ";

  parameter Modelica.Units.SI.ThermalResistance R_th=1/ka "K/W; Thermal resistance of stack";

  parameter Modelica.Units.SI.Voltage V_tn=1.482 "or U_tn, Thermoneutral voltage (voltage at which reaction can occur without releasing any heat)";

  parameter Modelica.Units.SI.Voltage v_n=simCenter.v_n "Nominal Voltage for grid";

  parameter Boolean usePowerPort=true "True if power port shall be used" annotation (Dialog(group="Replaceable Components"));

  parameter Boolean useHeatPort=true "True if heat port shall be used for cooling" annotation (Dialog(group="Replaceable Components"));

   //Temperature PID Controller Parameters
  parameter Modelica.Units.SI.Time tau_i=0.1 "1/tau_i for cooling system PID integrator gain";
  parameter Real k_p=980 "gain, cooling system PID proportional control - 1050 when opposite sign convention with PID";
  parameter Modelica.Units.SI.Time tau_d=1e-1 "tau_d, for cooling system PID derivator gain";
  parameter Real N_i=0.5 "gain of anti-windup compensation ";
  parameter Real N_d=1 "gain, ideal derivative block ";

  // _____________________________________________
  //
  //             Variable Declarations
  // _____________________________________________

  Modelica.Units.SI.Temperature T_cell=T_stack "Temperature of one cell";
  Modelica.Units.SI.Temperature T_stack(start=23.5 + 273.15) "Temperature of the stack" annotation (Dialog(group="Initialization", showStartAttribute=true));
  Modelica.Units.SI.Temperature T_syng_ein "Temperature of the syngas";
  Modelica.Units.SI.Temperature T_air_ein "Temperature of the air";

  Modelica.Units.SI.Pressure P_O2=air.p_i[3]/air.p "Partial pressure of oxygen at the cathode";
  Modelica.Units.SI.Pressure P_H2=syng.p_i[5]/syng.p "Partial pressure of hydrogen at the anode";

  Modelica.Units.SI.Voltage E_cell "Voltage of one cell";
  Modelica.Units.SI.Voltage E_stack "Voltage of the stack";
  Modelica.Units.SI.Voltage V_ohmic "Ohmic losses";
  Modelica.Units.SI.Voltage V_Nernst "Nernst voltage";
  Modelica.Units.SI.Voltage V_act "Activation voltage";
  //   Modelica.Units.SI.Voltage V_a "Dynamic activation voltage - Xue 2004 model";
  Modelica.Units.SI.Voltage V_rev "Reversible electric potential";

  Modelica.Units.SI.Conductivity mem_conductivity "Conductivity of the PE membrane";

  Real da= cp_tab[3,1] - cp_tab[1,1] - 0.5* cp_tab[2,1] "Empiric parameter for calculating Gibbs free energy according to Barbir";
  Real db= cp_tab[3,2] - cp_tab[1,2] - 0.5* cp_tab[2,2] "Empiric parameter for calculating Gibbs free energy according to Barbir";
  Real dc= cp_tab[3,3] - cp_tab[1,3] - 0.5* cp_tab[2,3] "Empiric parameter for calculating Gibbs free energy according to Barbir";

  Real Delta_H_T = -241.98*1000 + da*(T_cell-298.15) + db * ((T_cell^2) - 298.15^2)/2 + dc * ((T_cell^3) - 298.15^3)/3 "Empiric equation for calculating the enthalpy of formation according to Barbir";
  Real Delta_S_T = -0.0444*1000 + da*log(T_cell/298.15) + db * (T_cell - 298.15) + dc * ((T_cell^2) - 298.15^2)/2 "Empiric equation for calculating the entropy according to Barbir";

  Modelica.Units.SI.MolarFlowRate N_dot_e "Molar flow rate of the electrons";
  Modelica.Units.SI.MassFlowRate m_dot_H2_react_stack "Required H2 mass flow rate of one cell";
  Modelica.Units.SI.MassFlowRate m_dot_O2_react_stack "Required O2 mass flow rate of one cell";
  Modelica.Units.SI.MassFlowRate m_dot_H2O_gen_stack "Generated H2O mass flow rate of one cell";
  Modelica.Units.SI.MassFlowRate m_dot_air_react_stack "Required air mass flow rate of one cell";
  Modelica.Units.SI.MolarMass M_H2=syng.M_i[5] "Molar mass H2";
  Modelica.Units.SI.MolarMass M_O2=air.M_i[3] "Molar mass O2";
  Modelica.Units.SI.MassFraction xi_O2 "Mass fraction of O2 in the air";

  Modelica.Units.SI.SpecificEnthalpy h_hein=syng.h;
  Modelica.Units.SI.SpecificEnthalpy h_haus=synga.h;
  Modelica.Units.SI.SpecificEnthalpy h_aein=air.h;
  Modelica.Units.SI.SpecificEnthalpy h_aaus=aira.h;

  Boolean cooling_control "control operation of Q_cool";

  Modelica.Units.SI.HeatFlowRate Q_flow_reac "Heat flow due to reaction";
  Modelica.Units.SI.HeatFlowRate Q_flow_gas "Heat flow to/from syngas and air";
  Modelica.Units.SI.HeatFlowRate Q_flow_convective "Convective heat flow ";
  Modelica.Units.SI.HeatFlowRate Q_flow_cooling "cooling power of heat exchanger";
  Modelica.Units.SI.HeatFlowRate Q_flow_el "W, heat generated in FC from hydrogen reaction for electricity production";


  Modelica.Units.SI.Current I(start=1) "Electric current through the stack";
  Modelica.Units.SI.Current I_is "Theoretical??? electric current through the stack";
  Modelica.Units.SI.CurrentDensity i_cell "Current density";
  Modelica.Units.SI.Resistance Ri "Internal resistance in Ohm - intially Ri = - 1.667e-4*T_stack + 0.2289";
  // If use of dynamic activation overvoltage
   // Modelica.Units.SI.Resistance Ra "Activation resistance in Ohm";


  // for heating up states:
  Boolean is_Shutdown(start=false) "true, if load current is indicating to shut the stack down";

  // _____________________________________________
  //
  //                  Interfaces
  // _____________________________________________

  TransiEnt.Basics.Interfaces.Gas.IdealGasTempPortIn feedh(Medium=Syngas) annotation (Placement(transformation(extent={{-110,30},{-90,50}}), iconTransformation(extent={{-114,46},{-86,74}})));
  TransiEnt.Basics.Interfaces.Gas.IdealGasTempPortOut drainh(Medium=Syngas) annotation (Placement(transformation(extent={{90,30},{110,50}}), iconTransformation(extent={{86,46},{114,74}})));
  TransiEnt.Basics.Interfaces.Gas.IdealGasTempPortIn feeda(Medium=Air) annotation (Placement(transformation(extent={{-106,-66},{-86,-46}}), iconTransformation(extent={{-114,-74},{-86,-46}})));
  TransiEnt.Basics.Interfaces.Gas.IdealGasTempPortOut draina(Medium=Air) annotation (Placement(transformation(extent={{90,-64},{110,-44}}), iconTransformation(extent={{86,-74},{114,-46}})));

  replaceable TransiEnt.Basics.Interfaces.Electrical.ActivePowerPort epp if usePowerPort constrainedby TransiEnt.Basics.Interfaces.Electrical.PartialPowerPort "Choice of power port" annotation (Dialog(group="Replaceable Components"),choicesAllMatching=true, Placement(transformation(extent={{-10,48},{10,68}}), iconTransformation(extent={{-10,48},{10,68}})));

  TransiEnt.Basics.Interfaces.Electrical.ElectricCurrentIn I_load "Input for loading current" annotation (Placement(
        transformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={-50,100}),
                         iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=0,
        origin={-82,-6})));
  TransiEnt.Basics.Interfaces.Electrical.VoltageOut V_stack "Output for Voltage of a stack"  annotation (Placement(
        transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={50,100}),iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=0,
        origin={100,0})));
  TransiEnt.Basics.Interfaces.General.MassFractionOut lambda_H "Output for excess ratio of hydrogen" annotation (Placement(
        transformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={42,-100}),
                         iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={52,-100})));

  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_el = -1 * V_stack * I annotation (Placement(
        transformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={-6,-100}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=0,
        origin={106,-96})));
  TransiEnt.Basics.Interfaces.General.MassFractionOut lambda_O "Output for excess ratio of oxygen" annotation (Placement(
        transformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={-56,-102}),iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={-48,-100})));

  // _____________________________________________
  //
  //           Instances of other Classes
  // _____________________________________________

  TILMedia.Gas_pT syng(
    p=p_Anode,
    T=T_syng_ein,
    xi= inStream(feedh.xi_outflow),
    gasType = Syngas)
    annotation (Placement(transformation(extent={{-80,28},
            {-60,48}})));

    TILMedia.Gas_pT air(
    T=T_air_ein,
    p=p_Kathode,
    xi=inStream(feeda.xi_outflow),
    gasType = Air)
    annotation (Placement(transformation(extent={{62,24},{86,54}})));

    TILMedia.Gas_pT synga(
    p=p_Amb,
    T=T_stack,
    xi= drainh.xi_outflow,
    gasType = Syngas)
    annotation (Placement(transformation(extent={{-80,-52},
            {-60,-32}})));

    TILMedia.Gas_pT aira(
    T=T_stack,
    p=p_Amb,
    xi=draina.xi_outflow,
    gasType = Air)
    annotation (Placement(transformation(extent={{64,-72},{88,-42}})));

  replaceable TransiEnt.Components.Boundaries.Electrical.ActivePower.Power powerBoundary
                                                                               if usePowerPort constrainedby TransiEnt.Components.Boundaries.Electrical.Base.PartialModelPowerBoundary "Choice of power boundary model. The power boundary model must match the power port." annotation (
    Dialog(group="Replaceable Components"),
    choices(
      choice(redeclare TransiEnt.Components.Boundaries.Electrical.ActivePower.Power powerBoundary "P-Boundary for ActivePowerPort"),
      choice(redeclare TransiEnt.Components.Boundaries.Electrical.ApparentPower.ApparentPower powerBoundary(
          useInputConnectorP=true,
          useInputConnectorQ=false,
          useCosPhi=true,
          cosphi_boundary=1) "PQ-Boundary for ApparentPowerPort"),
      choice(redeclare TransiEnt.Components.Boundaries.Electrical.ComplexPower powerBoundary(useInputConnectorQ=false, cosphi_boundary=1) "PQ-Boundary for ComplexPowerPort"),
      choice(redeclare TransiEnt.Components.Boundaries.Electrical.ApparentPower.PowerVoltage powerBoundary(Use_input_connector_v=false, v_boundary=PEM.v_n) "PV-Boundary for ApparentPowerPort"),
      choice(redeclare TransiEnt.Components.Electrical.QuasiStationaryComponentsBusses.PUMachine powerBoundary(v_gen=PEM.v_n, useInputConnectorP=true) "PV-Boundary for ComplexPowerPort")),
    Placement(transformation(extent={{-22,18},{-42,38}})));
  Modelica.Blocks.Sources.RealExpression powerInput(y=P_el) if usePowerPort annotation (Placement(transformation(extent={{-62,54},{-42,74}})));


  Modelica.Blocks.Sources.RealExpression T_op_out(y=T_stack)
                                                          annotation (Placement(transformation(extent={{-44,-50},{-24,-30}})));
  TransiEnt.Basics.Interfaces.General.TemperatureOut temperatureOut annotation (Placement(transformation(extent={{-54,-40},{-34,-20}}), iconTransformation(extent={{-54,-40},{-34,-20}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow(T_ref=T_nom, alpha=0)
                                                                              if useHeatPort annotation (Placement(transformation(extent={{70,-36},{78,-28}})));
  Modelica.Blocks.Sources.RealExpression Q_flow_positive(y=-Q_flow_cooling)
                                                                     if useHeatPort annotation (Placement(transformation(extent={{50,-38},{62,-26}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heat if useHeatPort annotation (Placement(transformation(extent={{90,-42},{112,-20}}), iconTransformation(extent={{90,-42},{112,-20}})));
  Modelica.Blocks.Sources.RealExpression PID_T_max(y=T_stack_max)
                                                               annotation (Placement(transformation(extent={{-10,-48},{10,-28}})));
  Modelica.Blocks.Sources.BooleanExpression PID_T_control(y=cooling_control) annotation (Placement(transformation(extent={{-10,-62},{10,-42}})));
  Modelica.Blocks.Sources.RealExpression PID_T_op(y=T_stack)
                                                          annotation (Placement(transformation(extent={{-10,-76},{10,-56}})));
  ClaRa.Components.Utilities.Blocks.LimPID cooling_PID(
    y_max=0,
    y_min=-6000,
    Ni=N_i,
    y_inactive=0,
    use_activateInput=true,
    sign=1,
    Tau_i=tau_i,
    k=k_p,
    controllerType=Modelica.Blocks.Types.SimpleController.PI)  annotation (Placement(transformation(extent={{20,-56},{40,-36}})));
  Modelica.Blocks.Sources.RealExpression Q_Cooling(y=Q_flow_cooling) "excess waste heat generated by fuel cell system, actively cooled by default" annotation (Placement(transformation(extent={{-42,-80},{-22,-60}})));
  TransiEnt.Basics.Interfaces.Thermal.HeatFlowRateOut excessHeatFlowOut annotation (Placement(transformation(extent={{-54,-70},{-34,-50}})));

equation
  // _____________________________________________
  //
  //           Characteristic equations
  // _____________________________________________

  // Current and voltage equations
  // Equations from Fuel cell technology by N. Sammes, p. 228
  // Equations from System level lumped-parameter dynamic modeling of PEM fuel cell, X. Xue (2004) and Amphlett (1996)
  // Activation voltage coefficients come from Amphlett (1996), A model predicting transient responses of PEMFC

  V_rev = 1.229 - 8.5e-4*(T_stack - T_amb);
  V_Nernst = ( Modelica.Constants.R * T_stack) / (z * Modelica.Constants.F) * log(P_H2 * sqrt(P_O2));
  V_ohmic = - I*t_mem / (mem_conductivity * A_cell);
  //   V_ohmic = - I * Ri; // results in a negative voltage at high power setpoints
  mem_conductivity = (0.00514*lambda - 0.00326)*exp(1268*(1/298.15 - 1/T_stack))*10e2; // Springer, 1991 "Polymer Electrolyte Fuel Cell Model"
  Ri = 0.016046 - 3.4715e-5*T_stack + 7.9565e-5*I;

  E_cell = V_rev + V_Nernst + V_act + V_ohmic;
  E_stack = E_cell * no_Cells;

// If use of dynamic activation overvoltage
//   E_cell = V_rev + V_Nernst + V_a + V_ohmic;
//   der(V_a) = 2 * I - 2 * V_a * Ra;
//   Ra = - V_act * I;


  I_is = 2*feedh.m_flow*inStream(feedh.xi_outflow[5])*Modelica.Constants.F / (M_H2*no_Cells);
  i_cell = I / A_cell;

  // Reactions mass flows
  N_dot_e = I / Modelica.Constants.F;

  m_dot_H2_react_stack =  N_dot_e / 2 * M_H2 * no_Cells;
  m_dot_O2_react_stack =  N_dot_e / 4 * M_O2 * no_Cells;
  m_dot_H2O_gen_stack = N_dot_e / 2 * syng.M_i[4] * no_Cells;
  xi_O2 = 1 - (air.xi[1] + air.xi[2]);
  m_dot_air_react_stack = m_dot_O2_react_stack / xi_O2;


 // Shutting down model
  is_Shutdown = I_load < I_shutdown;

  if is_Shutdown then
    // Load current below mimimum value = signal to shut down!
    I = 0;
    V_stack = 54;  // in Xue (2004), V_cell tends towards 36.75 V for I = 0. Based on results from these simulations, highest reached voltage at 5000 W (OCV approximation) is closer to 54 V. We will only see the OCV when P_el = 0 and I = 0. Allowed power setpoints are either 0 or greater than 500 Wm and computed by MPC
    V_act = 0;
    lambda_H = -1; // signal to lambda h controller that we are in shut down mode!
    lambda_O = -1;
  else
    // Normal operating point: Reaction is running
    //     I = min(I_load,I_is);
    I = I_load;
    V_stack = E_stack;
    //     V_act = -0.9514 + 0.00312*T_stack + 7.4e-5*T_stack*log(I) - 1.87e-4*T_stack* log(P_O2*1.97e5*exp(398.15/T_stack)); // Amphlett 1995
    V_act = -0.9514 + 0.00312*T_stack + 7.4e-5*T_stack*log(I) - 1.87e-4*T_stack* log(xi_O2); // xi_O2 is a mass fraction, but interior of log () must be in mol/cm3 - we multiply xi_O2 by O2 density (0.001225 g/cm3) divided by O2 molecular weight (32 g/mol) - less coherent results with that
    lambda_H = (feedh.m_flow*inStream(feedh.xi_outflow[5]))/m_dot_H2_react_stack;
    lambda_O = (feeda.m_flow*(1-inStream(feeda.xi_outflow[1])-inStream(feeda.xi_outflow[2])))/m_dot_O2_react_stack;
  end if;



  // temperature design flow direction
  T_syng_ein = inStream(feedh.T_outflow);
  T_air_ein = inStream(feeda.T_outflow);
  drainh.T_outflow = synga.T;
  draina.T_outflow = aira.T;

  // contra design flow direction
  feeda.T_outflow = aira.T;
  feedh.T_outflow = aira.T;

  // Energy balance
  cooling_control = T_stack > T_cool_set;

  Q_flow_el = if (E_cell < V_tn) then 0 else no_Cells*(E_cell - V_tn)*I;

  Q_flow_reac = -1*m_dot_H2O_gen_stack/syng.M_i[4]*Delta_H_T;

  Q_flow_gas = feedh.m_flow*h_hein + feeda.m_flow*h_aein + drainh.m_flow*h_haus + draina.m_flow*h_aaus;

  Q_flow_convective = (T_amb - T_stack)/R_th "convection heat transfer rate";

  Q_flow_cooling = cooling_PID.y;


  der(T_stack)*m*cp = Q_flow_reac + Q_flow_gas + Q_flow_convective + Q_flow_cooling - Q_flow_el;



  // mass balance (total mass)
  feeda.m_flow - m_dot_air_react_stack = - draina.m_flow;
  feedh.m_flow - m_dot_H2_react_stack = - drainh.m_flow;

  // mass fractions (design flow direction)
  -1*draina.m_flow*draina.xi_outflow[1] = feeda.m_flow*inStream(feeda.xi_outflow[1]) + m_dot_H2O_gen_stack/no_Cells;  //Water mass flow rate of cell stack because of the inserted mass and the reaction
  -1*draina.m_flow*draina.xi_outflow[2] = feeda.m_flow*inStream(feeda.xi_outflow[2]);  //Nitrogen mass flow rate of the cell stack - if balanced for H2O and N, balanced for O2
  -1*drainh.m_flow*drainh.xi_outflow[1] = feedh.m_flow*inStream(feedh.xi_outflow[1]);
  -1*drainh.m_flow*drainh.xi_outflow[2] = feedh.m_flow*inStream(feedh.xi_outflow[2]);
  -1*drainh.m_flow*drainh.xi_outflow[3] = feedh.m_flow*inStream(feedh.xi_outflow[3]);
  -1*drainh.m_flow*drainh.xi_outflow[4] = feedh.m_flow*inStream(feedh.xi_outflow[4]);
  -1*drainh.m_flow*drainh.xi_outflow[5] = feedh.m_flow*inStream(feedh.xi_outflow[5]) - m_dot_H2_react_stack/no_Cells;   //Hydrogen mass flow rate of the cell stack
  -1*drainh.m_flow*drainh.xi_outflow[6] = feedh.m_flow*inStream(feedh.xi_outflow[6]);

  // mass fraction (opposite design flow direction)
  feeda.xi_outflow[1] = 0;
  feeda.xi_outflow[2] = 0;
  for i in 1:6 loop
  feedh.xi_outflow[i] = 0;
  end for;

  // impulse equation
  feedh.p = drainh.p;
  feeda.p = draina.p;


  if usePowerPort then
  connect(powerBoundary.epp, epp) annotation (Line(
      points={{-22,28},{-22,27.95},{0,27.95},{0,58}},
      color={0,127,0},
      thickness=0.5));
    connect(powerInput.y, powerBoundary.P_el_set) annotation (Line(points={{-41,64},{-28,64},{-28,40},{-26,40}}, color={0,0,127}));
  end if;

  if useHeatPort then
  connect(heat,prescribedHeatFlow.port) annotation (Line(points={{101,-31},{90,-31},{90,-32},{78,-32}},
                                                                                            color={191,0,0}));
  end if;

  connect(temperatureOut,T_op_out. y) annotation (Line(points={{-44,-30},{-18,-30},{-18,-40},{-23,-40}},
                                                                                    color={0,0,127}));
  connect(PID_T_op.y, cooling_PID.u_m) annotation (Line(points={{11,-66},{30.1,-66},{30.1,-58}},                  color={0,0,127}));
  connect(PID_T_control.y, cooling_PID.activateInput) annotation (Line(points={{11,-52},{12,-54},{18,-54}},      color={255,0,255}));
  connect(PID_T_max.y, cooling_PID.u_s) annotation (Line(points={{11,-38},{14,-38},{14,-46},{18,-46}},
                                                                                               color={0,0,127}));
  connect(excessHeatFlowOut, Q_Cooling.y) annotation (Line(
      points={{-44,-60},{-18,-60},{-18,-70},{-21,-70}},
      color={162,29,33},
      pattern=LinePattern.Dash));

  connect(Q_flow_positive.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{62.6,-32},{70,-32}}, color={0,0,127}));
  connect(heat, heat) annotation (Line(points={{101,-31},{101,-31}}, color={191,0,0}));
 annotation (Line(points={{50.6,-26},{70,-26}}, color={0,0,127}),
              Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics={
        Rectangle(
          extent={{-54,44},{-34,-48}},
          lineColor={95,95,95},
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{30,44},{50,-48}},
          lineColor={95,95,95},
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-34,44},{-22,-48}},
          lineColor={95,95,95},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{18,44},{30,-48}},
          lineColor={95,95,95},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-22,44},{-10,-48}},
          lineColor={95,95,95},
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{6,44},{18,-48}},
          lineColor={95,95,95},
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-10,44},{6,-48}},
          lineColor={95,95,95},
          fillColor={255,0,128},
          fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
                Documentation(info="<html>
<h4><span style=\"color: #008000\">1. Purpose of model</span></h4>
<p>Model of PEM cell stack.</p>
<h4><span style=\"color: #008000\">2. Level of detail, physical effects considered, and physical insight</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">3. Limits of validity </span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">4. Interfaces</span></h4>
<p>Type of electrical power port can be chosen</p>
<h4><span style=\"color: #008000\">5. Nomenclature</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">6. Governing Equations</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">7. Remarks for Usage</span></h4>
<p><span style=\"font-family: MS Shell Dlg 2;\">With the choice of the boundary the the model can be used as PQ or PU bus.</span></p>
<h4><span style=\"color: #008000\">8. Validation</span></h4>
<p><span style=\"font-family: MS Shell Dlg 2;\">Validation has been done as part of the master thesis [1] </span></p>
<p><b><span style=\"font-family: MS Shell Dlg 2; color: #008000;\">9. References</span></b></p>
<p><span style=\"font-family: MS Shell Dlg 2;\">[1] </span>Modellierung und Simulation von erdgasbetriebenen Brennstoffzellen-Blockheizkraftwerken zur Heimenergieversorgung</p>
<p>Master thesis, Simon Weilbach (2014) </p>
<h4><span style=\"color: #008000\">10. Version History</span></h4>
<p>Model created by Simon Weilbach (simon.weilbach@tuhh.de) <span style=\"font-family: MS Shell Dlg 2;\">in October 2014</span></p>
<p><span style=\"font-family: MS Shell Dlg 2;\">Model revised by Pascal Dubucq (dubucq@tuhh.de) in October 2015</span></p>
<p><span style=\"font-family: MS Shell Dlg 2;\">Quality check (Code conventions) by Rebekka Denninger in October 2016</span></p>
<p><span style=\"font-family: MS Shell Dlg 2;\">Model generalized for different electrical power ports by Jan-Peter Heckel (jan.heckel@tuhh.de) in July 2018 </span></p>
</html>"));
end PEMFC;
