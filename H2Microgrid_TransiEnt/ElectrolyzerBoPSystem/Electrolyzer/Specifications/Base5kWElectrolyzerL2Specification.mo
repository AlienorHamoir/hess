within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.Specifications;
record Base5kWElectrolyzerL2Specification "Record used for specification of an Electrolyzer system"

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

  extends TransiEnt.Basics.Icons.Record;

  import      Modelica.Units.SI;

  // _____________________________________________
  //
  //                   Parameters
  // _____________________________________________

  parameter Integer n_cells=20 "number of cells in series in the electrolyzer stack" annotation(Dialog(group="Fundamental Definitions"));
  parameter Real alpha_an=0.8 "charge transfer coefficient for anode" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.Energy E_exc=52994 "Activation energy for anode reaction" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.ChemicalPotential E_pro=10542 "Temperature independent parameter for activation energy in proton transport" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.Area PEM_area=50e-4 "Active surface area of electrolyzer cell electrodes" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.Current i_el_n=75 "Nominal current of electrolyzer" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.ActivePower P_el_n=5.5e3 "Nominal power of the electrolyzer" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.ActivePower P_el_max=2.0*P_el_n "Maximum power of the electrolyzer" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.ThermalResistance R_th=0.529056 "Thermal resistance of stack" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.HeatCapacity C_th=54038.66667 "Overall lumped thermal capacitance" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.CurrentDensity i_dens_0_an_std=1e-7 "Exchange current density at electrodes at T_std" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.Thickness t_mem=254e-6 "PE membrane thickness" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.Conductivity mem_conductivity_ref=10.47 "Membrane conductivity value at T_std" annotation(Dialog(group="Fundamental Definitions"));
  parameter Real specificWaterConsumption=10 "Mass of water per mass of hydrogen" annotation(Dialog(group="Fundamental Definitions")); //Stolzenburg, K. et al.: Integration von Wind-Wasserstoff-Systemen in das Energiesystem: Abschlussbericht, 2014
  parameter SI.Temperature T_out=273.15 + 40 "Hydrogen output temperature" annotation(Dialog(group="Fundamental Definitions"));
  parameter SI.Power P_el_pump = 285.7 "pump el power consumption" annotation(Dialog(group="Cooling circuit"));
  parameter SI.Efficiency eta_pumpmotor=0.63 "pump's motor electric efficiency" annotation(Dialog(group="Cooling circuit"));
  parameter SI.VolumeFlowRate V_flow_water= 0.000127
                                                    "water flow rate in cooling loop" annotation(Dialog(group="Cooling circuit"));
  parameter SI.PressureDifference Delta_p_pump=6.9*100000 "total pump head" annotation(Dialog(group="Cooling circuits"));
  parameter SI.Temperature T_op_max=273.15 + 75 "max operating temp" annotation(Dialog(group="Cooling circuit"));
  parameter SI.Temperature T_cool_set=273.15 + 50 "Cooling trigger point" annotation(Dialog(group="Cooling circuit"));
  parameter SI.HeatFlowRate Q_flow_cool_max=6911 "Maximum cooling power" annotation(Dialog(group="Cooling circuit"));

  annotation (Documentation(info="<html>
<h4><span style=\"color: #008000\">1. Purpose of model</span></h4>
<p>Partial record for electrolyzer system specific data.</p>
<h4><span style=\"color: #008000\">2. Level of detail, physical effects considered, and physical insight</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">3. Limits of validity </span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">4. Interfaces</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">5. Nomenclature</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">6. Governing Equations</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">7. Remarks for Usage</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">8. Validation</span></h4>
<p>(no validation or testing necessary)</p>
<h4><span style=\"color: #008000\">9. References</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">10. Version History</span></h4>
<p>Adopted from BaseCHPRecord. Created by John Webster (jcwebste@edu.uwaterloo.ca) Oct. 2018</p>
</html>"));
end Base5kWElectrolyzerL2Specification;
