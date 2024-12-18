within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer;
model PEMElectrolyzer_L2 "PEMElectrolyzer_L2 Proton exchange membrane electrolyzer with selectable physics submodels, including inverter efficiency, compressor and cooling pump power, total power setpoint, system and stack efficiencies"

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

 extends TransiEnt.Producer.Gas.Electrolyzer.Base.PartialElectrolyzer(P_el_n=Specification.P_el_n,P_el_max=Specification.P_el_max,T_out=Specification.T_out,
    realExpression(y=P_el_tot),
    Q_flow_positive(y=-Q_flow_heatprovision),
    prescribedHeatFlow(T_ref=T_op_n, alpha=0.0001),
    realExpression8(y=T_op_n),
    heatFlow_externalMassFlowControl(
      change_sign=false,
      T_out_limit_const=348.15,
      useVariableToutlimit=true));

  // _____________________________________________
  //
  //        Constants and Hidden Parameters
  // _____________________________________________

  inner parameter Modelica.Units.SI.Temperature T_std=298.15 "STD temperature";

  // _____________________________________________
  //
  //              Visible Parameters
  // _____________________________________________

  //Temperature parameters
public
  inner parameter Modelica.Units.SI.Temperature T_amb=296.65        "K, ambient temperature" annotation (Dialog(group="Fundamental Definitions"));
  parameter Modelica.Units.SI.Temperature T_op_start=273.15 + 23.5 "initial operating temperature of PEM Electrolyzer" annotation (Dialog(group="Initialization"));
  parameter Modelica.Units.SI.Temperature T_op_n=273.15 + 50 "nominal operating temperature of PEM Electrolyzer";

  //Electrolyzer system specific parameters
protected
  inner parameter Integer n_cells=Specification.n_cells "Number of electrolysis cells in series";
  inner parameter Real alpha_an=Specification.alpha_an "Charge transfer coefficient for anode";
  inner parameter Modelica.Units.SI.Energy E_exc=Specification.E_exc "Activation energy for anode reaction";
  inner parameter Modelica.Units.SI.ChemicalPotential E_pro=Specification.E_pro "Temperature independent parameter for activation energy in proton transport";
  inner parameter Modelica.Units.SI.CurrentDensity i_dens_0_an_std=Specification.i_dens_0_an_std "Exchange current density at electrodes at T_std";
  inner parameter Modelica.Units.SI.Thickness t_mem=Specification.t_mem "PE membrane thickness";
  inner parameter Modelica.Units.SI.Conductivity mem_conductivity_ref=Specification.mem_conductivity_ref "Membrane conductivity value at T_std";
  inner parameter Modelica.Units.SI.Area PEM_area=Specification.PEM_area "Active surface area of electrolyzer cell electrodes";
  inner parameter Modelica.Units.SI.Current i_el_n=Specification.i_el_n "Nominal current of electrolyzer";
  inner parameter Modelica.Units.SI.ThermalResistance R_th=Specification.R_th "Thermal resistance of stack";
  inner parameter Modelica.Units.SI.HeatCapacity C_th=Specification.C_th "Lumped thermal capacitance of stack";
  parameter Real specificWaterConsumption=Specification.specificWaterConsumption "Mass of water per mass of hydrogen";
  inner parameter Modelica.Units.SI.Power P_el_pump=Specification.P_el_pump "pump el power consumption when electrolyzer is running";
  inner parameter Modelica.Units.SI.Efficiency eta_pumpmotor=Specification.eta_pumpmotor "pump's motor electric efficiency";
  inner parameter Modelica.Units.SI.VolumeFlowRate V_flow_water=Specification.V_flow_water "water flow rate in cooling loop";
  inner parameter Modelica.Units.SI.PressureDifference Delta_p_pump=Specification.Delta_p_pump "total pump head";
  inner parameter Modelica.Units.SI.Temperature T_op_max=Specification.T_op_max "Max operating temp";
  inner parameter Modelica.Units.SI.Temperature T_cool_set=Specification.T_cool_set "Cooling trigger point";
  inner parameter Modelica.Units.SI.HeatFlowRate Q_flow_cool_max=Specification.Q_flow_cool_max "Maximum cooling power";

  //Model Configuration Parameters
public
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3 "Medium model" annotation (Dialog(group="Fundamental Definitions"));
  inner parameter Integer whichInput=3 "Use current, current density, electric power, or mass_flow_H2 as input" annotation (Dialog(group="Fundamental Definitions"), choices(
      __Dymola_radioButtons=true,
      choice=1 "Current",
      choice=2 "Current Density",
      choice=3 "Electric Power",
      choice=4 "mass_flow_H2"));
  inner parameter Boolean userSetTemp=false "Use T_input as input" annotation (Dialog(group="Fundamental Definitions"), choices(
      __Dymola_radioButtons=true,
      choice=true "yes",
      choice=false "no"));
//  parameter Real eta_inv_n=0.956 "nominal efficiency of the inverter" annotation (Dialog(group="Fundamental Definitions")); // we neglect all losses relative to grid constraints, including inverter efficiency
  parameter Real E_dry_spec=1400*3600
                                    "specific energy consumption of dryer in Ws/kg H2" annotation (Dialog(group="Fundamental Definitions"));

  // _____________________________________________
  //
  //              Variable Declarations
  // _____________________________________________

public
  Modelica.Units.SI.Mass mass_H2(start=0, fixed=false) "Produced H2 mass";

  //Temperature modeling
  //The following must all be calculated in the Temperature model or else provided externally.
  inner Modelica.Units.SI.Temperature T_op "Operating stack temperature";

  //Power and energy, efficiency
  Modelica.Units.SI.Power P_el "Electric power consumed by the electrolyzer";
  Modelica.Units.SI.Power P_circ_pump "Electric power consumed by the circulation pump";
  Modelica.Units.SI.Power P_el_tot "Electric power consumed by the electrolyzer and dryer and cooling and water circulation";
  Modelica.Units.SI.Power P_el_sys "Electric power actually consumed by the system electrolyzer and dryer and cooling and water circulation";
  Modelica.Units.SI.Power P_el_aux "Electric power consumed by the electrolyzer auxiliaries system";


  Modelica.Units.SI.Energy E_dry "Electric energy consumed by the dryer";
  //   Modelica.Units.SI.Efficiency eta_cond "Neglected here as already included in drying energy - Variable for modeling the efficiency loss of the electrolyseur due to the needed power for the dryer and the water conditioning";
  Modelica.Units.SI.EnthalpyFlowRate H_flow_H2_GCV "H2 enthalpy flow rate out of electrolyzer, gross calorific value";
  Modelica.Units.SI.EnthalpyFlowRate H_flow_H2_NCV "H2 enthalpy flow rate out of electrolyzer, net calorific value";
  Modelica.Units.SI.Efficiency eta_NCV_EL(max=1, min=0) "Efficiency of the electrolyzer based on NCV" annotation (Dialog(group="Initialization", showStartAttribute=true));
  Modelica.Units.SI.Efficiency eta_GCV_EL(max=1, min=0) "Efficiency of the electrolyzer based on GCV" annotation (Dialog(group="Initialization", showStartAttribute=true));
  Modelica.Units.SI.Efficiency eta_NCV_sys(max=1, min=0) "Efficiency of the electrolyzer + BoP system based on NCV" annotation (Dialog(group="Initialization", showStartAttribute=true));
  Modelica.Units.SI.Efficiency eta_GCV_sys(max=1, min=0) "Efficiency of the electrolyzer + BoP system based on GCV" annotation (Dialog(group="Initialization", showStartAttribute=true));

  //Current
  inner Modelica.Units.SI.Current i_el_stack "Current across the electrolyzer stack";
  inner Modelica.Units.SI.CurrentDensity i_dens_a(start=0.01, fixed=false) "Operating current density at anode";

  //Voltage and Overpotential Variables
  //The following must all be calculated in the Voltage model.
  inner Modelica.Units.SI.Voltage V_el_stack "PEM stack voltage" annotation (Dialog(group="Fundamental Definitions"));
  inner Modelica.Units.SI.Voltage V_cell "PEM cell voltage considering all included physical phenomena";
  inner Modelica.Units.SI.Voltage V_tn "or U_tn, Thermoneutral voltage (voltage at which reaction can occur without releasing any heat)";

  //Pressure modeling
  //The following must all be calculated in the Pressure model
  inner Modelica.Units.SI.Pressure gasPortPressure "pressure of hydrogen connected to gasPortOut";
  inner Modelica.Units.SI.Pressure pp_H2O "Pa, vapour pressure of water vapour, must always be converted to atm";
  inner Modelica.Units.SI.Pressure pp_H2 "Pa, partial pressure of H2, must always be converted to atm";
  inner Modelica.Units.SI.Pressure pp_O2 "Pa, partial pressure of O2, must always be converted to atm";
  inner Modelica.Units.SI.Pressure p_cat "Pressure at cathode (H2)";
  inner Modelica.Units.SI.Pressure p_an "Pressure at anode (O2)";
  //   inner SI.Pressure  p_sat "water saturation pressure";

  //Mass Transport modeling
  //The following must all be calculated in the Mass Flow model, or provided externally (m_flow_H2).
  inner Modelica.Units.SI.MassFlowRate m_flow_H2 "H2 mass flow rate out of electrolyzer";
  inner Modelica.Units.SI.MolarFlowRate n_flow_H2O "Molar consummation rate of water";
  inner Modelica.Units.SI.MolarFlowRate n_flow_H2 "Molar production rate of Hydrogen gas";
  inner Modelica.Units.SI.MolarFlowRate n_flow_O2 "Molar production rate of Oxygen gas";

  model Outline
    extends TransiEnt.Basics.Icons.Record;
    input Modelica.Units.SI.Voltage V_cell "PEME cell voltage";
    input Modelica.Units.SI.Voltage V_stack "PEME stack voltage";
    input Modelica.Units.SI.Current i_stack "PEME current consumption";
    input Modelica.Units.SI.Temperature T_op "ambient temperature of electrolyzer operation";
    input Modelica.Units.SI.HeatFlowRate der_Q_cool "heat to coolant system";
    input Modelica.Units.SI.HeatFlowRate der_Q_loss "heat loss in PEM stack";
    input Modelica.Units.SI.HeatFlowRate W_pemwe "heat generated in PEM stack";
    input Modelica.Units.SI.Power P_el "Consumed electric power";
    input Modelica.Units.SI.Energy W_el "Consumed electric energy";
    input Modelica.Units.SI.Power H_flow_NCV "Produced enthalpy flow based on NCV";
    input Modelica.Units.SI.Power H_flow_GCV "Produced enthalpy flow based on GCV";
    input Modelica.Units.SI.Mass mass_H2 "Produced hydrogen mass";
    input Modelica.Units.SI.Mass mass_H2O "Consumed water mass";
    input Modelica.Units.SI.Efficiency eta_NCV "Efficiency based on NCV";
    input Modelica.Units.SI.Efficiency eta_GCV "Efficiency based on GCV";
    /*
    input SI.HeatFlowRate W_pump_BHP;
    input SI.HeatFlowRate W_pump_elec;
    input SI.HeatFlowRate W_pump_WHP;
    input SI.HeatFlowRate W_pump;
    input SI.HeatFlowRate W_lost;
    */
    /*
    input SI.Voltage V_ohmic "test ohmic OP";
    input SI.Voltage V_activation "test electrode OP";
    input SI.Voltage V_rev "rev cell voltage";
    input SI.Voltage V_nernst "nernst voltage";
    input Real R_mem "Membrane resistance";
    input Real mem_conductivity "mem_conductivity";
    input Real i_a "current density on anode";
    input Real i_0a "exchange current density on anode";
    */
  end Outline;

  model Summary
    extends TransiEnt.Basics.Icons.Record;
    Outline outline;
    TransiEnt.Basics.Records.FlangeRealGas gasPortOut;
    TransiEnt.Basics.Records.Costs costs;
  end Summary;

  // _____________________________________________
  //
  //                 Outer Models
  // _____________________________________________

  replaceable model electrolyzerVoltage = TransiEnt.Producer.Gas.Electrolyzer.Base.Physics.Voltage.V_cell1
                                                                                                      constrainedby TransiEnt.Producer.Gas.Electrolyzer.Base.Physics.Voltage.PartialV_cell
                                                                                                                                                                                      "Dynamic voltage behaviour of electrolyzer" annotation (
    Dialog(group="Fundamental Definitions"),
    choicesAllMatching=true,
    Placement(transformation(extent={{-30,-10},{-10,10}})));
  replaceable model electrolyzerTemperature = TransiEnt.Producer.Gas.Electrolyzer.Base.Physics.Temperature.Temperature1
                                                                                                                   constrainedby TransiEnt.Producer.Gas.Electrolyzer.Base.Physics.Temperature.PartialTemperature
                                                                                                                                                                                                        "Temperature model" annotation (
    Dialog(group="Fundamental Definitions"),
    Placement(transformation(extent={{10,-10},{30,10}})),
    choicesAllMatching=true);

  replaceable model electrolyzerPressures = TransiEnt.Producer.Gas.Electrolyzer.Base.Physics.Pressures.Pressures1
                                                                                                             constrainedby TransiEnt.Producer.Gas.Electrolyzer.Base.Physics.Pressures.PartialPressures
                                                                                                                                                                                                  "Pressure model" annotation (
    Dialog(group="Fundamental Definitions"),
    Placement(transformation(extent={{50,-10},{10,10}})),
    choicesAllMatching=true);

  replaceable model electrolyzerMassFlow = TransiEnt.Producer.Gas.Electrolyzer.Base.Physics.MassFlow.MassFlow0thOrderDynamics
                                                                                               constrainedby TransiEnt.Producer.Gas.Electrolyzer.Base.Physics.MassFlow.PartialMassFlow
                                                                                                                                                                                  "Mass flow model" annotation (
    Dialog(group="Fundamental Definitions"),
    Placement(transformation(extent={{90,-10},{50,10}})),
    choicesAllMatching=true);

  // _____________________________________________
  //
  //           Instances of other Classes
  // _____________________________________________

  replaceable H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.Specifications.GinerELX5kW Specification constrainedby H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.Specifications.Base5kWElectrolyzerL2Specification "Record containing technical data of electrolyzer system" annotation (
    Placement(transformation(extent={{-14,-98},{6,-78}})),
    Dialog(tab="General", group="Specification"),
    choices(choicesAllMatching=true));

  electrolyzerVoltage voltage annotation (Placement(transformation(extent={{-22,-2},{-8,12}})));

  electrolyzerTemperature temperature(cooling_PID(controllerType=Modelica.Blocks.Types.SimpleController.PI,
      y_max=Q_flow_cool_max))         annotation (Placement(transformation(extent={{-20,-26},{-6,-12}})));

  electrolyzerPressures pressure annotation (Placement(transformation(extent={{16,-2},{2,12}})));

  electrolyzerMassFlow massFlow annotation (Placement(transformation(extent={{16,-26},{2,-12}})));

  // _____________________________________________
  //
  //                Interfaces
  // _____________________________________________

  Modelica.Blocks.Interfaces.RealInput i_el_stack_set(displayUnit="A", final unit="A")
                if whichInput==1 annotation (Placement(transformation(
        extent={{-14,-14},{14,14}},
        rotation=270,
        origin={22,94}), iconTransformation(
        extent={{-14,-14},{14,14}},
        rotation=270,
        origin={22,94})));

  Modelica.Blocks.Interfaces.RealInput i_dens_set(displayUnit="A/m2", final unit="A/m2")
                       if whichInput==2 annotation (Placement(transformation(
        extent={{-28,-14},{-1.77636e-15,14}},
        rotation=270,
        origin={2,80}), iconTransformation(
        extent={{-28,-14},{-1.77636e-15,14}},
        rotation=270,
        origin={2,80})));

  Modelica.Blocks.Interfaces.RealInput P_el_set(
    final quantity="Power",
    displayUnit="W",
    final unit="W") if whichInput==3 annotation (Placement(transformation(
        extent={{-14,-14},{14,14}},
        rotation=270,
        origin={-46,94}), iconTransformation(
        extent={{-14,-14},{14,14}},
        rotation=270,
        origin={-46,94})));

  Modelica.Blocks.Interfaces.RealInput m_flow_H2_set(
    final quantity="MassFlowRate",
    displayUnit="kg/s",
    final unit="kg/s") if whichInput==4 annotation (Placement(transformation(
        extent={{-14,-14},{14,14}},
        rotation=270,
        origin={46,94}), iconTransformation(
        extent={{-14,-14},{14,14}},
        rotation=270,
        origin={46,94})));

  Modelica.Blocks.Interfaces.RealInput T_input(
    final quantity="Temperature",
    displayUnit="K",
    final unit="K") if userSetTemp annotation (Placement(transformation(
        extent={{-14,-14},{14,14}},
        rotation=270,
        origin={-22,94}), iconTransformation(
        extent={{-14,-14},{14,14}},
        rotation=270,
        origin={-22,94})));

  TransiEnt.Basics.Interfaces.Gas.RealGasPortOut gasPortOut(Medium=medium) annotation (Placement(transformation(extent={{90,-10},{110,10}})));

protected
  Modelica.Blocks.Sources.Constant i_el_stack_set_(k=0) if not whichInput==1;
  Modelica.Blocks.Sources.Constant i_dens_set_(k=0) if not whichInput==2;
  Modelica.Blocks.Sources.Constant P_el_set_(k=0) if not whichInput==3;
  Modelica.Blocks.Sources.Constant m_flow_H2_set_(k=0) if not whichInput==4;

  Modelica.Blocks.Sources.Constant T_input_(k=0) if not userSetTemp;

  TransiEnt.Producer.Gas.Electrolyzer.Base.GetInputsElectrolyzerL2 getInputs annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=270,
        origin={0,58})));

public
  inner Summary summary(
    outline(
      V_cell=V_cell,
      V_stack=V_el_stack,
      i_stack=i_el_stack,
      T_op=T_op,
      der_Q_cool=temperature.Q_flow_cooling,
      der_Q_loss=temperature.Q_flow_convective,
      W_pemwe=temperature.Q_flow_H2Oelectrolysis,
      P_el=P_el,
      W_el=collectCosts.W_el_demand,
      mass_H2=mass_H2,
      mass_H2O=specificWaterConsumption*mass_H2,
      eta_NCV=eta_NCV_sys,
      eta_GCV=eta_GCV_sys,
      H_flow_NCV=gasPortOut.m_flow*NCV_H2[end],
      H_flow_GCV=gasPortOut.m_flow*GCV_H2[end]),
    gasPortOut(
      mediumModel=medium,
      xi=vleFluidH2.xi,
      x=vleFluidH2.x,
      m_flow=-gasPortOut.m_flow,
      T=vleFluidH2.T,
      p=gasPortOut.p,
      h=vleFluidH2.h,
      rho=vleFluidH2.d),
    costs(
      costs=collectCosts.costsCollector.Costs,
      investCosts=collectCosts.costsCollector.InvestCosts,
      demandCosts=collectCosts.costsCollector.DemandCosts,
      oMCosts=collectCosts.costsCollector.OMCosts,
      otherCosts=collectCosts.costsCollector.OtherCosts,
      revenues=collectCosts.costsCollector.Revenues)) annotation (Placement(transformation(extent={{10,-98},{30,-78}})));

  /* //additional outline variables
      W_pump_BHP=temperature.W_pump_BHP,
      W_lost=temperature.W_lost,
      W_pump_elec=temperature.W_pump_elec,
      W_pump_WHP=temperature.W_pump_WHP,
      W_pump=temperature.W_pump,

      i_a=i_dens_a,
      i_0a=voltage.i_dens_0a,
      mem_conductivity=voltage.mem_conductivity,
      R_mem=voltage.R_mem,
      V_ohmic=voltage.V_ohmic,
      V_activation=voltage.V_activation,

      V_rev=voltage.V_rev,
      V_nernst=voltage.V_nernst
*/

  TransiEnt.Basics.Interfaces.General.TemperatureOut temperatureOut "Stack operating temperature"
                                                                    annotation (Placement(transformation(extent={{-82,-38},{-62,-18}}), iconTransformation(extent={{-82,-38},{-62,-18}})));
  Modelica.Blocks.Sources.RealExpression T_op_out(y=T_op) annotation (Placement(transformation(extent={{-100,-34},{-80,-14}})));
  TransiEnt.Basics.Interfaces.Thermal.HeatFlowRateOut excessHeatFlowOut "Excess heat flow that has to be cooled to maintain adequate stack temperature "
                                                                        annotation (Placement(transformation(extent={{-82,-58},{-62,-38}}), iconTransformation(extent={{-82,-58},{-62,-38}})));
  Modelica.Blocks.Sources.RealExpression Q_Flow_Cooling(y=temperature.Q_flow_cooling) "excess waste heat generated by electrolyzer system, actively cooled by default" annotation (Placement(transformation(extent={{-100,-52},{-80,-32}})));

  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn CoolingPumpPowerIn "Electrical power from the cooling system pump"
                                                                            annotation (Placement(transformation(
        extent={{-14,-16},{14,16}},
        rotation=180,
        origin={114,16}), iconTransformation(
        extent={{-14,-16},{14,16}},
        rotation=180,
        origin={114,16})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn CompressorPowerIn "Electrical power from the storage compressor" annotation (Placement(transformation(
        extent={{-13,-15},{13,15}},
        rotation=180,
        origin={113,43}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=0,
        origin={-120,-76})));
  Modelica.Blocks.Sources.RealExpression P_el_out(y=P_el_sys) annotation (Placement(transformation(extent={{-100,-68},{-80,-48}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut electrolyzerPowerOut "Net power consumed by the electrolyzer and its auxiliary systems"
                                                                               annotation (Placement(transformation(extent={{-82,-80},{-62,-60}}), iconTransformation(extent={{-82,-80},{-62,-60}})));
  Modelica.Blocks.Sources.RealExpression P_el_out_aux(y=P_el_aux) annotation (Placement(transformation(extent={{-100,-84},{-80,-64}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut electrolyzerPowerOutAux "Net power consumed by the electrolyzer auxiliary systems" annotation (Placement(transformation(extent={{-82,-100},{-62,-80}}), iconTransformation(extent={{-82,-100},{-62,-80}})));
initial equation
  if userSetTemp == false then
    T_op = T_op_start;
  end if;
equation
  // _____________________________________________
  //
  //           Characteristic Equations
  // ___________________________________________
  //GasPortOut
  gasPortOut.xi_outflow = xi_out;
  gasPortOut.h_outflow = vleFluidH2.h;
  gasPortOut.m_flow = -m_flow_H2;

  gasPortPressure = gasPortOut.p;
  Q_flow=-(temperature.Q_flow_cooling);

  if whichInput == 1 then
    //set current
    i_el_stack = getInputs.i_el_stack_set;
  elseif whichInput == 2 then
    //set current density
    i_dens_a = getInputs.i_dens_set;
  elseif whichInput == 3 then
    //set electric power
    P_el_tot = getInputs.P_el_set; // eventually, P_el_setpoint should take into account consumption of whole electrolyzer system
    //P_el = getInputs.P_el_set; // used for system calibration
  elseif whichInput == 4 then
    //set mass flow
    m_flow_H2 = getInputs.m_flow_H2_set;
  end if;

  if userSetTemp then
    T_op = getInputs.T_input;
  end if;

  if integrateH2Flow then
    der(mass_H2) = m_flow_H2;
  else
    mass_H2 = 0;
  end if;

  m_flow_H2O = specificWaterConsumption*m_flow_H2;

  //Power
  V_el_stack = V_cell*n_cells;
  i_el_stack = i_dens_a*PEM_area;
  P_el = (V_el_stack)*i_el_stack;
  E_dry = E_dry_spec*mass_H2; // drying energy = 1.4 kWh/kg H2
  P_el_tot =P_el + der(E_dry) + CoolingPumpPowerIn + CompressorPowerIn + P_circ_pump;
  P_el_aux =der(E_dry) + CoolingPumpPowerIn + CompressorPowerIn + P_circ_pump;


  if P_el>0 then
    P_el_sys = P_el_tot;
    P_circ_pump = P_el_pump;
  else
    P_el_sys = 0;
    P_circ_pump = 0;
  end if;

  //Efficiency

  // Neglect the eta_cond term as drying (and water conditioning is neglected) is already taken into account in E_dry
  H_flow_H2_GCV = (m_flow_H2*(GCV_H2[end] + (vleFluidH2.h - h0)));
  H_flow_H2_NCV = (m_flow_H2*(NCV_H2[end] + (vleFluidH2.h - h0)));

  if P_el > 0 then
    eta_GCV_EL = H_flow_H2_GCV  /P_el;
    eta_NCV_EL = H_flow_H2_NCV  /P_el;
    eta_GCV_sys = H_flow_H2_GCV /P_el_tot;
    eta_NCV_sys = H_flow_H2_NCV  /P_el_tot;
  else
    eta_GCV_EL = 0;
    eta_NCV_EL = eta_GCV_EL;
    eta_GCV_sys = 0;
    eta_NCV_sys = eta_GCV_sys;
  end if;

  // _____________________________________________
  //
  //           Connect Statements
  // _____________________________________________

  connect(i_el_stack_set_.y, getInputs.i_el_stack_set);
  connect(i_dens_set_.y, getInputs.i_dens_set);
  connect(P_el_set_.y, getInputs.P_el_set);
  connect(m_flow_H2_set_.y, getInputs.m_flow_H2_set);

  connect(T_input_.y, getInputs.T_input);

  connect(i_el_stack_set, getInputs.i_el_stack_set) annotation (Line(points={{22,94},{22,74},{3.6,74},{3.6,70}},       color={0,0,127}));
  connect(i_dens_set, getInputs.i_dens_set) annotation (Line(points={{2,94},{1,92},{1,70}},                  color={0,0,127}));
  connect(m_flow_H2_set, getInputs.m_flow_H2_set) annotation (Line(points={{46,94},{46,70},{6.6,70}},                color={0,0,127}));

  connect(T_input, getInputs.T_input) annotation (Line(points={{-22,94},{-22,74},{-1.4,74},{-1.4,70}},  color={0,0,127}));

  connect(temperatureOut, T_op_out.y) annotation (Line(points={{-72,-28},{-76,-28},{-76,-24},{-79,-24}},
                                                                                    color={0,0,127}));
  connect(excessHeatFlowOut, Q_Flow_Cooling.y) annotation (Line(
      points={{-72,-48},{-68,-48},{-68,-42},{-79,-42}},
      color={162,29,33},
      pattern=LinePattern.Dash));
  connect(P_el_set, getInputs.P_el_set) annotation (Line(points={{-46,94},{-46,70},{-4,70}},          color={0,0,127}));
  connect(P_el_out.y, electrolyzerPowerOut) annotation (Line(points={{-79,-58},{-68,-58},{-68,-70},{-72,-70}},
                                                                                           color={0,0,127}));
  connect(excessHeatFlowOut, excessHeatFlowOut) annotation (Line(
      points={{-72,-48},{-72,-48}},
      color={175,0,0},
      pattern=LinePattern.Dash));
  connect(P_el_out_aux.y, electrolyzerPowerOutAux) annotation (Line(points={{-79,-74},{-68,-74},{-68,-90},{-72,-90}}, color={0,0,127}));
  annotation (defaultComponentName="electrolyzer", Documentation(info="<html>
<h4><span style=\"color: #008000\">1. Purpose of model</span></h4>
<p>This model electrolyzer uses modular physics classes and a Specification record to describe a real-life electrolyzer system. The default model uses physical relationships taken from (Espinosa-L&oacute;pez et al, 2018).</p>
<h4><span style=\"color: #008000\">2. Level of detail, physical effects considered, and physical insight</span></h4>
<p>The user can select or create a model with system specific parameters. The desired input (electric power, current, current density or hydrogen mass flow) can be varied. Physics submodels can be replaced as desired, with essential definitions annotated. The water consumption, hydrogen output, and net- and gross-calorific energy conversion efficiency can be calculated. The user has the option of controlling pressure through gasPortOut and/or temperature through T_input as well.</p>
<p>It has been modified from inital TransiEnt model :</p>
<p>(1) so it is able to account for storage compressor power and cooling pump power from the cooling system if desired. The electrical power setpoint accounts for the entire electrolyzer system consumption, including auxiliaries consumption.</p>
<p>(2) temperature PID has been tuned against 5kW Giner electrolyzer parameters</p>
<p>(3) initial term for water conditioning power has been neglected; the oly conditionig is the dryig of the hydroge, which is considered proportional to the mass flow rate of produced hydrogen</p>
<h4><span style=\"color: #008000\">3. Limits of validity </span></h4>
<p>Original model developed and validated in the range of 20-60 &deg;C with operating pressure of 15-35 bar. </p>
<h4><span style=\"color: #008000\">4. Interfaces</span></h4>
<p>epp: electric power port, type can be chosen </p>
<p>gasPortOut: hydrogen outlet </p>
<p>i_dens_set: input for electric current density</p>
<p>i_el_stack_set: input for electric current</p>
<p>P_el_set: input for electric power </p>
<p>m_flow_H2_set: input for hydrogen mass flow </p>
<p>T_input: input for temperature</p>
<p>excessHeatFlowOut: Heat flow rate out port, equal to the cooling power used to regulate temperature at max temp.</p>
<p>temperatureOut: temperatureOut interface equal to operating temperature</p>
<h4><span style=\"color: #008000\">5. Nomenclature</span></h4>
<p>(no elements)</p>
<h4><span style=\"color: #008000\">6. Governing Equations</span></h4>
<p>Selectable physics equations allow for different governing equations to be used, and consist of equations from (Espinosa-L&oacute;pez et al, 2018) by default.</p>
<h4><span style=\"color: #008000\">7. Remarks for Usage</span></h4>
<p>(no remarks) </p>
<h4><span style=\"color: #008000\">8. Validation</span></h4>
<p>Results have been validated against (Espinosa-L&oacute;pez et al, 2018) published figures and Abdin results. </p>
<h4><span style=\"color: #008000\">9. References</span></h4>
<p>Source: modified from TransiEnt library original model</p>
<p><br>[1] Manuel Espinosa-L&oacute;pez, Philippe Baucour, Serge Besse, Christophe Darras, Raynal Glises, Philippe Poggi, Andr&eacute; Rakotondrainibe, and Pierre Serre-Combe. Modelling and experimental validation of a 46 kW PEM high pressure water electrolyser. Renewable Energy, 119, pp. 160-173, 2018. doi: 10.1016/J.RENENE.2017.11.081. </p>
<p>[2] efficiency curve of the inverter taken from the data sheet of SMA &quot;Technische Wirkungsgrade und Derating&quot; URL: https://files.sma.de/dl/1348/WirkungDerat-TI-de-46.pdf page 71, 26.11.2019</p>
<p>[3] J. Webster and C. Bode, &ldquo;Implementation of a Non-Discretized Multiphysics PEM Electrolyzer Model in Modelica,&rdquo; in Proceedings of the 13th International Modelica Conference, Regensburg, Germany, March 4&ndash;6, 2019, no. 157, pp. 833&ndash;840, DOI: 10.3384/ecp19157833.</p>
<p>[4] Z. Abdin, E. MacA. Gray, and C.J. Webb. Modelling and simulation of a proton exchange membrane (PEM) electrolyzer cell. International Journal of Hydrogen Energy, 40(39):13243-13257, 2015. doi:10.1016/j.ijhydene.2015.07.129. </p>
<h4><span style=\"color: #008000\">10. Version History</span></h4>
<p>Model created by John Webster (jcwebste@edu.uwaterloo.ca) in October 2018</p>
<p>Model adjusted for TransiEnt by Jan Westphal (j.westphal@tuhh.de) in dec 2019</p>
</html>"));
end PEMElectrolyzer_L2;
