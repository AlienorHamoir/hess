within H2Microgrid_TransiEnt.FuelCellBoPSystem.FuelCell;
model PEMFC_KhanIqbal "Model of PEM Fuel Cell stack experimentally validated - Khan and Iqbal"

//________________________________________________________________________________//
// Adapted from PEMFC component of the TransiEnt Library, version: 2.0.3                             //
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

  parameter Integer no_Cells = 36 "Number of cells connected in series";

  parameter Real lambda=11.1074 "constant humidity";

  parameter Modelica.Units.SI.Thickness t_mem = 125e-6 "PE membrane thickness - ref. Dow's membrane in microm";

  parameter Real z = 2 "Quantity of transfered electrons";

  parameter Modelica.Units.SI.Area A=0.0232 "Area of one cell";

  parameter Modelica.Units.SI.Pressure p_Anode=3.03975e5 "Pressure at the anode - 3 atm";

  parameter Modelica.Units.SI.Pressure p_Cathode=3.03975e5  "Pressure at the cathode - 3 atm";

  parameter Modelica.Units.SI.Pressure p_Amb=1e5 "Ambient pressure - 1 bar";

  parameter Modelica.Units.SI.Pressure p_Atm=1.01325e5 "Atmospheric pressure - 1 atm";

  parameter TransiEnt.Basics.Media.Gases.Gas_VDIWA_SG7_var Syngas=TransiEnt.Basics.Media.Gases.Gas_VDIWA_SG7_var() "Medium model of Syngas" annotation (choicesAllMatching);

  parameter TransiEnt.Basics.Media.Gases.Gas_MoistAir Air=TransiEnt.Basics.Media.Gases.Gas_MoistAir() "Medium model of air - Moist air is used which consists of N2, H2O and O2. O2 is used in FC" annotation (choicesAllMatching);

  parameter Real cp_tab[3,3] = [28.91404, -0.00084, 2.01e-6; 25.84512, 0.012987, -3.9e-6; 30.62644, 0.009621, 1.18e-6] "Empiric parameter for calculating the Gibbs free energy according to Barbir";

  parameter Modelica.Units.SI.Mass m=43 "Mass of the stack";

  parameter Modelica.Units.SI.SpecificHeatCapacity cp=814 "Specific heat capacity of the stack - cp * m = 35000";

  parameter Modelica.Units.SI.Temperature T_nom = 72 + 273.15 "Temperature in nominal point";

  parameter Modelica.Units.SI.Temperature T_amb = 23.5 + 273.15 "Ambient temperature";

  parameter Modelica.Units.SI.Temperature T_std = 25 + 273.15 "Standard temperature";

  parameter Modelica.Units.SI.Temperature T_stack_max = 75 + 273.15 "Maximum stack temperature";

  parameter Modelica.Units.SI.Temperature T_cool_set = 73 + 273.15 "Cooling trigger point";

  parameter Modelica.Units.SI.ThermalConductance ka=22 "Thermal conductance of the stack ";

  parameter Modelica.Units.SI.ThermalResistance R_th=1/ka "°C/W; Thermal resistance of stack";

  parameter Modelica.Units.SI.Voltage V_tn=1.482 "or U_tn, Thermoneutral voltage (voltage at which reaction can occur without releasing any heat)";

  parameter Modelica.Units.SI.Voltage v_n=simCenter.v_n "Nominal Voltage for grid";

  parameter Modelica.Units.SI.Voltage OCV=0.99 "Open circuit voltage of a cell";

  parameter Modelica.Units.SI.Voltage E_0=0.89195 "Open circuit voltage of a cell used in voltage equation";

  parameter Modelica.Units.SI.Capacitance C_dl = 4  "Stack double layer capacitance for dynamic activation overvoltage";

  parameter Modelica.Units.SI.Enthalpy H_0 = 285.5e3 "Hydrogen enthalpy of combustion";

  parameter Modelica.Units.SI.SpecificHeatCapacity cp_H2O = 4182 "Specific heat capacity of liquid water";

  parameter Modelica.Units.SI.Current I_shutdown=8 "If load controller requests currents below this value, stack will shut down";

  parameter Modelica.Units.SI.Volume V_a = 0.005 "Anode volume";

  parameter Modelica.Units.SI.Volume V_c = 0.01 "Cathode volume";

  parameter Real k_a = 0.065 "Anode flow constant in mol/ s atm";

  parameter Real k_c = 0.065 "Cathode flow constant in mol/ s atm";

  parameter Boolean usePowerPort=false "True if power port shall be used" annotation (Dialog(group="Replaceable Components"));

  parameter Boolean useHeatPort=false "True if heat port shall be used for cooling" annotation (Dialog(group="Replaceable Components"));

  parameter Real eps = 1e-6 "Tolerance to avoid mathematical errors";

   //Temperature PID Controller Parameters
  parameter Modelica.Units.SI.Time tau_i=0.01 "1/tau_i, cooling system PID integrator gain";
  parameter Real k_p=3800 "gain, cooling system PID proportional control";
  parameter Modelica.Units.SI.Time tau_d=0.1 "tau_d, for cooling system PID derivator gain";
  parameter Real N_i=0.5 "gain of anti-windup compensation ";
  parameter Real N_d=1 "gain, ideal derivative block ";

  // _____________________________________________
  //
  //             Variable Declarations
  // _____________________________________________
  Real homotopyFactor = homotopy(actual=1, simplified=0);

  Modelica.Units.SI.Temperature T_cell=T_stack "Temperature of one cell";
  Modelica.Units.SI.Temperature T_stack(start=T_amb, fixed=true) "Temperature of the stack" annotation (Dialog(group="Initialization", showStartAttribute=true));
  Modelica.Units.SI.Temperature T_syng_ein "Temperature of the syngas";
  Modelica.Units.SI.Temperature T_air_ein "Temperature of the air";

  Modelica.Units.SI.Voltage E_cell "Voltage of one cell";
  Modelica.Units.SI.Voltage E_stack "Voltage of the stack";
  Modelica.Units.SI.Voltage V_ohmic "Ohmic losses";
  Modelica.Units.SI.Voltage V_Nernst "Nernst voltage";
  Modelica.Units.SI.Voltage V_act "Activation voltage";
 // Modelica.Units.SI.Voltage V_d "Dynamic activation voltage - Correa dynamic model";
  Modelica.Units.SI.Voltage V_conc "Concentration voltage";
  Modelica.Units.SI.Resistance Rm "Membrane resistance in Ohm";
  // If use of dynamic activation overvoltage
  //Modelica.Units.SI.Resistance Ra "Activation resistance in Ohm";

  Real da= cp_tab[3,1] - cp_tab[1,1] - 0.5* cp_tab[2,1] "Empiric parameter for calculating Gibbs free energy according to Barbir";
  Real db= cp_tab[3,2] - cp_tab[1,2] - 0.5* cp_tab[2,2] "Empiric parameter for calculating Gibbs free energy according to Barbir";
  Real dc= cp_tab[3,3] - cp_tab[1,3] - 0.5* cp_tab[2,3] "Empiric parameter for calculating Gibbs free energy according to Barbir";

  Real Delta_H_T = -241.98*1000 + da*(T_cell-298.65) + db * ((T_cell^2) - 298.65^2)/2 + dc * ((T_cell^3) - (298.65^3))/3 "Empiric equation for calculating the enthalpy of formation according to Barbir";
  Real Delta_S_T = -0.0444*1000 + da*log(T_cell/298.15) + db * (T_cell - 298.15) + dc * ((T_cell^2) - 298.15^2)/2 "Empiric equation for calculating the entropy according to Barbir";

  Modelica.Units.SI.MolarFlowRate N_dot_e "Molar flow rate of the electrons";
  Modelica.Units.SI.MassFlowRate m_dot_H2_react_stack "Required H2 mass flow rate of one cell";
  Modelica.Units.SI.MassFlowRate m_dot_O2_react_stack "Required O2 mass flow rate of one cell";
  Modelica.Units.SI.MassFlowRate m_dot_H2O_gen_stack "Generated H2O mass flow rate of one cell";
  Modelica.Units.SI.MassFlowRate m_dot_air_react_stack "Required air mass flow rate of one cell";
  Modelica.Units.SI.MolarMass M_H2=syng.M_i[5] "Molar mass H2";
  Modelica.Units.SI.MolarMass M_H2O=syng.M_i[4] "Molar mass H2O";
  Modelica.Units.SI.MolarMass M_O2=air.M_i[3] "Molar mass O2";
  Modelica.Units.SI.MassFraction xi_O2 "Mass fraction of O2 in the air";
  Modelica.Units.SI.MoleFraction x_O2 "Molar fraction of O2 in the air";
  Modelica.Units.SI.MassFraction xi_H2 "Mass fraction of H2 in the air";
  Modelica.Units.SI.MoleFraction x_H2 "Molar fraction of H2 in the air";
  Modelica.Units.SI.Concentration c_O2 "Concentration in mol/m3 fraction of O2 in the stack";
  Modelica.Units.SI.Concentration c_H2 "Concentration in mol/m3 fraction of H2 in the stack";
  Modelica.Units.SI.Resistivity r_mem "Resistivity of the PE membrane";
  Modelica.Units.SI.Pressure P_O2=x_O2*air.p "Partial pressure of oxygen at the cathode";
  Modelica.Units.SI.Pressure P_H2=x_H2*syng.p "Partial pressure of hydrogen at the anode";
  // Modelica.Units.SI.Pressure P_H2O "Partial pressure of water at the cathode";

  Modelica.Units.SI.SpecificEnthalpy h_hein=syng.h;
  Modelica.Units.SI.SpecificEnthalpy h_haus=synga.h;
  Modelica.Units.SI.SpecificEnthalpy h_aein=air.h;
  Modelica.Units.SI.SpecificEnthalpy h_aaus=aira.h;

  Boolean cooling_control "control operation of Q_cool";
  Modelica.Units.SI.TemperatureSlope der_T_stack "Derivative of the operating stack temperature";
  Modelica.Units.SI.HeatFlowRate Q_flow_gas "Heat flow from syngas";
  Modelica.Units.SI.HeatFlowRate Q_flow_convective "Convective heat flow ";
  Modelica.Units.SI.HeatFlowRate Q_flow_cooling "Cooling power of heat exchanger";
  Modelica.Units.SI.HeatFlowRate Q_flow_el "W, heat used in FC by hydrogen reaction for electricity production";

  Boolean is_Shutdown(start=false) "true, if load current is indicating to shut the stack down";
  Modelica.Units.SI.Current I "Electric current through the stack";
  Modelica.Units.SI.Current I_is "Theoretical electric current through the stack based on available H2 mass flow";
  Modelica.Units.SI.CurrentDensity J "Current density";

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

  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_el annotation (Placement(
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
    p=p_Cathode,
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
  Modelica.Blocks.Sources.RealExpression Q_flow_positive(y=Q_flow_cooling)
                                                                     if useHeatPort annotation (Placement(transformation(extent={{50,-38},{62,-26}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heat if useHeatPort annotation (Placement(transformation(extent={{90,-42},{112,-20}}), iconTransformation(extent={{90,-42},{112,-20}})));
  Modelica.Blocks.Sources.BooleanExpression PID_T_control(y=cooling_control) annotation (Placement(transformation(extent={{-10,-62},{10,-42}})));
  Modelica.Blocks.Sources.RealExpression PID_T_op(y=T_stack)
                                                          annotation (Placement(transformation(extent={{-10,-76},{10,-56}})));
  ClaRa.Components.Utilities.Blocks.LimPID cooling_PID(
    y_max=15000,
    y_min=0,
    Tau_d=tau_d,
    Ni=N_i,
    Nd=N_d,
    y_inactive=0,
    use_activateInput=true,
    sign=-1,
    Tau_i=tau_i,
    k=k_p,
    controllerType=Modelica.Blocks.Types.SimpleController.PI)  annotation (Placement(transformation(extent={{20,-56},{40,-36}})));
  Modelica.Blocks.Sources.RealExpression Q_Cooling(y=Q_flow_cooling) "excess waste heat generated by fuel cell system, actively cooled by default" annotation (Placement(transformation(extent={{-42,-80},{-22,-60}})));
  TransiEnt.Basics.Interfaces.Thermal.HeatFlowRateOut excessHeatFlowOut annotation (Placement(transformation(extent={{-54,-70},{-34,-50}})));

  Modelica.Blocks.Sources.RealExpression T_setpoint(y=T_nom) annotation (Placement(transformation(extent={{-10,-44},{10,-24}})));
equation
  // _____________________________________________
  //
  //           Characteristic equations
  // _____________________________________________

          ////////
  //// Voltage and current model ////
          ////////
  P_el = - V_stack * I;

  V_Nernst = E_0 - 0.9e-3*(T_stack - T_amb) + ( Modelica.Constants.R * T_stack) / (z * Modelica.Constants.F) * log(P_H2 *1e-5 * sqrt(P_O2 * 1e-5)); // replace E0 by OCV = 0.99 V for more accurate polarization curve

  c_O2 = P_O2 * 1.97e-7 * exp(498/T_stack);
  c_H2 = P_H2 * 9.174e-7 * exp(-77/T_stack);

  V_ohmic = I * Rm;
  Rm =   t_mem * r_mem / A;
  r_mem = 181.6 / ((lambda - 0.634) * exp(4.18 * ((T_stack-303)/T_stack)))* 1e-2; // simplified Khan and Iqbal (2004), neglecting dependency on J

  V_conc = 0;

  // If use of steady-state activation overvoltage
  E_cell = V_Nernst - V_act - V_ohmic - V_conc;

  // If use of dynamic activation overvoltage
  //der(V_d) = I/C_dl - V_d / (C_dl*Ra);
  //Ra =  V_act / I;
  //E_cell = V_Nernst - V_d - V_ohmic - V_conc;

  E_stack = E_cell * no_Cells;

  I_is = 2*feedh.m_flow*inStream(feedh.xi_outflow[5])*Modelica.Constants.F / (M_H2*no_Cells);
  J = I / A; // A/m^2

 // Minimal load current management
  is_Shutdown = I_load < I_shutdown;

  if is_Shutdown then
    // Load current below mimimum value = signal to shut down/stand-by!
    I = 0;
    V_stack = no_Cells * 1.05; // replace E_0 by OCV
    V_act = 0;
    lambda_H = -1; // signal to lambda h controller that we are in shut down mode!
    lambda_O = -1;
  else
   // Normal operating point: Reaction is running
  //I = min(I_load,I_is); // to accomodate with the hydrogen supply dynamics
  I = I_load;
  V_stack = E_stack;
  V_act = 0.948 - (0.00286 + 0.0002 * log(A*1e4) + 4.3e-5 * log(c_H2*10e-6)) * T_stack - 7.6e-5*T_stack*log(c_O2*10e-6) + 1.93e-4*T_stack* log(I/(A*1e4));
  lambda_H = (feedh.m_flow*inStream(feedh.xi_outflow[5]))/m_dot_H2_react_stack;
  lambda_O = (feeda.m_flow*(1-inStream(feeda.xi_outflow[1]) - inStream(feeda.xi_outflow[2])))/m_dot_O2_react_stack;
  end if;

  N_dot_e = I / Modelica.Constants.F;

          ////////
  //// Thermal model ////
          ////////

 // temperature design flow direction
  T_syng_ein = inStream(feedh.T_outflow);
  T_air_ein = inStream(feeda.T_outflow);
  drainh.T_outflow = synga.T;
  draina.T_outflow = aira.T;

  feeda.T_outflow = air.T;
  feedh.T_outflow = syng.T;

  // Energy balance
  der(T_stack) = der_T_stack;

  cooling_control = T_stack > T_cool_set;

  Q_flow_gas = m_dot_H2_react_stack / M_H2 * H_0 + feedh.m_flow*h_hein + feeda.m_flow*h_aein + drainh.m_flow*h_haus + draina.m_flow*h_aaus - m_dot_H2O_gen_stack/M_H2O*cp_H2O*(T_stack-T_amb);

  Q_flow_el = E_stack*I;

  Q_flow_convective = (T_stack - T_amb)/R_th "convection heat transfer rate";

  Q_flow_cooling = cooling_PID.y;

  m*cp * der_T_stack = Q_flow_gas  - Q_flow_el - Q_flow_convective - Q_flow_cooling;

          ////////
  //// Reactant flow model ////
          ////////

   //// H2 balance at anode (syngas balance)

  // Dynamic model
//   V_a * der(P_H2) / (Modelica.Constants.R * T_stack) = (feedh.m_flow*inStream(feedh.xi_outflow[5]) - m_dot_H2_react_stack + drainh.m_flow*drainh.xi_outflow[5])/M_H2;
//   P_H2 = p_Anode - P_H2O;
//   - drainh.m_flow*drainh.xi_outflow[5] = M_H2 * k_a * (P_H2 - p_Anode);

  // Static model
  feedh.m_flow*inStream(feedh.xi_outflow[5]) = m_dot_H2_react_stack - drainh.m_flow*drainh.xi_outflow[5]; // H2 balance
  feedh.m_flow - m_dot_H2_react_stack = - drainh.m_flow; // syngas feed

  -drainh.m_flow*drainh.xi_outflow[1] = feedh.m_flow*inStream(feedh.xi_outflow[1]); // Methane balance
  -drainh.m_flow*drainh.xi_outflow[2] = feedh.m_flow*inStream(feedh.xi_outflow[2]); // Oxygen balance
  -drainh.m_flow*drainh.xi_outflow[3] = feedh.m_flow*inStream(feedh.xi_outflow[3]); // Carbon dioxide balance
  -drainh.m_flow*drainh.xi_outflow[4] = feedh.m_flow*inStream(feedh.xi_outflow[4]); // Water balance at anode - we assume all the water is formed at the cathode
  -drainh.m_flow*drainh.xi_outflow[6] = feedh.m_flow*inStream(feedh.xi_outflow[6]); // Carbon monoxide balance

  m_dot_H2_react_stack =  N_dot_e / 2 * M_H2 * no_Cells;
  xi_H2 = syng.xi[5];
  x_H2 = syng.x[5];

  //// O2 balance at cathode

  // Dynamic model
//   V_c * der(P_O2) / (Modelica.Constants.R * T_stack)  = (feeda.m_flow*(1-inStream(feeda.xi_outflow[1])-inStream(feeda.xi_outflow[2])) - m_dot_O2_react_stack + draina.m_flow*(1-draina.xi_outflow[1]-draina.xi_outflow[2]))/M_O2;
//   P_O2 = p_Cathode - P_H2O;
//   - draina.m_flow*(1-draina.xi_outflow[1]-draina.xi_outflow[2]) = M_O2 * k_c * (P_O2 - p_Cathode);

  // Static model
  -draina.m_flow*draina.xi_outflow[1] = feeda.m_flow*inStream(feeda.xi_outflow[1]) + m_dot_H2O_gen_stack;  //Water balance at cathode
  -draina.m_flow*draina.xi_outflow[2] = feeda.m_flow*inStream(feeda.xi_outflow[2]);  //Nitrogen mass flow rate of the cell stack - if balanced for H2O and N, balanced for O2
  feeda.m_flow - m_dot_air_react_stack = - draina.m_flow; // air feed

  m_dot_O2_react_stack =  N_dot_e / 4 * M_O2 * no_Cells;
  xi_O2 = 1 - (air.xi[1] + air.xi[2]);
  x_O2 = 1 - (air.x[1] + air.x[2]);
  m_dot_air_react_stack = m_dot_O2_react_stack / xi_O2;

  m_dot_H2O_gen_stack = N_dot_e / 2 * M_H2O * no_Cells;

  //// Products balance: H2O mass balance at cathode (since water is produced at cathode side)

  // Dynamic model
//   V_c * der(P_H2O) / (Modelica.Constants.R * T_stack)  = (feeda.m_flow*inStream(feeda.xi_outflow[1]) + m_dot_H2O_gen_stack + draina.m_flow*draina.xi_outflow[1])/ M_H2O;
//   P_H2O = exp(-2.1794 + 0.02953 * T_stack - 9.1837e-5 * T_stack^2 + 1.4454e-7*T_stack^7);
//   - draina.m_flow*draina.xi_outflow[1] = M_H2O * k_c * (P_H2O - p_Cathode);

  // mass fractions (design flow direction)
  //-draina.m_flow*draina.xi_outflow[1] - drainh.m_flow*drainh.xi_outflow[4] = feeda.m_flow*inStream(feeda.xi_outflow[1]) + feedh.m_flow*inStream(feedh.xi_outflow[4]) + m_dot_H2O_gen_stack;  //Water balance of the stack because of the inserted mass and the reaction

  // mass fraction (opposite design flow direction)
  feeda.xi_outflow[1] = 0;
  feeda.xi_outflow[2] = 0;
  for i in 1:6 loop
  feedh.xi_outflow[i] = 0;
  end for;

  // impulse equation
  feedh.p = p_Anode;
  feeda.p = p_Cathode;

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
  connect(excessHeatFlowOut, Q_Cooling.y) annotation (Line(
      points={{-44,-60},{-18,-60},{-18,-70},{-21,-70}},
      color={162,29,33},
      pattern=LinePattern.Dash));

  connect(Q_flow_positive.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{62.6,-32},{70,-32}}, color={0,0,127}));
  connect(heat, heat) annotation (Line(points={{101,-31},{101,-31}}, color={191,0,0}));
  connect(T_setpoint.y, cooling_PID.u_s) annotation (Line(points={{11,-34},{16,-34},{16,-38},{18,-38},{18,-46}}, color={0,0,127}));
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
<p>Model of low temperature PEM fuel cell stack.</p>
<h4><span style=\"color: #008000\">2. Level of detail, physical effects considered, and physical insight</span></h4>
<p>Dynamic temperature and mass flows models. </p>
<p>Quasi-steady electrochemical model: voltage dependency on temperature, with Nernst, ohmic, reversible and activation voltages contribution.</p>
<p>Steady pressure model.</p>
<p>From initial model, following modifications have been made:</p>
<ul>
<li>high temperature fuel cell adpated to a low temperature fuel cell: adapted eergy balance and temperature model, cooling PID added to allow connection with a cooling system (based on electrolyzer L2 model &quot;H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.PEMElectrolyzer_L2&quot;)</li>
<li>voltage equations have been adapted to better coincides with litterature and possibility of calibration with experimental parameters</li>
<li>syngas composition has been adapted (in PEMFC system) to correspond to pure H2</li>
<li>air and hydrogen mass balance has been corrected: to allow for a varying air intake and adapted syngas-to-H2 composition</li>
</ul>
<h4><span style=\"color: #008000\">3. Limits of validity </span></h4>
<p>Model parameters are adapted to a Ballard Mark V 35-cell 5 kW PEM fuel cell stack [1] [2] [3], and fine-tuned with simulation results.</p>
<p>A distinction is made between shutdow and normal operating point modes.</p>
<p>Temperature range of 20-80 &deg;C.</p>
<h4><span style=\"color: #008000\">4. Interfaces</span></h4>
<p>Type of electrical power port can be chosen</p>
<h4><span style=\"color: #008000\">5. Nomenclature</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">6. Governing Equations</span></h4>
<p>current and voltage equations</p>
<p>reaction mass flows equations </p>
<p>temperature design flow direction and contra design flow direction</p>
<p>energy balance</p>
<p>mass balance and mass fractions (design and opposite flow direction) </p>
<p>impulse equations (pressure conservation)</p>
<h4><span style=\"color: #008000\">7. Remarks for Usage</span></h4>
<p><span style=\"font-family: MS Shell Dlg 2;\">With the choice of the boundary the model can be used as PQ or PU bus.</span></p>
<h4><span style=\"color: #008000\">8. Validation</span></h4>
<p><span style=\"font-family: MS Shell Dlg 2;\">Initial model validation has been done as part of the master thesis [4] </span></p>
<p><span style=\"font-family: MS Shell Dlg 2;\">No experimental validation has been done on the adapted model.</span></p>
<p><b><span style=\"font-family: MS Shell Dlg 2; color: #008000;\">9. References</span></b></p>
<p>Source: initial model PEMFC from TransiEnt library &quot;TransiEnt.Components.Electrical.FuelCellSystems.FuelCell.PEM&quot;</p>
<p>[1] X. Xue, J. Tang, A. Smirnova et al. System level lumped-parameter dyamic modeling of PEM fuel cell. Journal of Power Sources, 133(2): 188-204. doi: 10.1016/j.jpowsour.2003.12.064.</p>
<p>[2] J.C. Amphlett, R.F. Mann, B.A. Peppley, P.R. Roberge, A. Rodrigues. A model predicting transient responses of proton exchange membrane fuel cells. (1996) Journal of Power Sources, 61 (1-2): 183-188. doi: 10.1016/S0378-7753(96)02360-9. </p>
<p>[3] Barbir. Fuel cell technology: Reaching towards commercialization. Engineering Materials and Processes.</p>
<p>[4] Master thesis, Simon Weilbach (2014). Modellierung und Simulation von erdgasbetriebenen Brennstoffzellen-Blockheizkraftwerken zur Heimenergieversorgung</p>
<h4><span style=\"color: #008000\">10. Version History</span></h4>
<p>Model created by Simon Weilbach (simon.weilbach@tuhh.de) <span style=\"font-family: MS Shell Dlg 2;\">in October 2014</span></p>
<p><span style=\"font-family: MS Shell Dlg 2;\">Model revised by Pascal Dubucq (dubucq@tuhh.de) in October 2015</span></p>
<p><span style=\"font-family: MS Shell Dlg 2;\">Quality check (Code conventions) by Rebekka Denninger in October 2016</span></p>
<p><span style=\"font-family: MS Shell Dlg 2;\">Model generalized for different electrical power ports by Jan-Peter Heckel (jan.heckel@tuhh.de) in July 2018 </span></p>
<p><span style=\"font-family: MS Shell Dlg 2;\">Model adapted for a low-temperature PEM fuel cell by Ali&eacute;nor Hamoir in June 2024</span></p>
</html>"),
    experiment(
      StartTime=1,
      StopTime=1000,
      Interval=1,
      Tolerance=1e-06,
      __Dymola_Algorithm="Rkfix4"));
end PEMFC_KhanIqbal;
