within H2Microgrid_TransiEnt;
model SimCenter "SimCenter for global parameters, ambient conditions and collecting statistics"

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

  extends ClaRa.SimCenter;
  extends TransiEnt.Basics.Icons.Grids;

  // _____________________________________________
  //
  //                   Components
  // _____________________________________________

  inner replaceable TransiEnt.Components.Boundaries.Ambient.AmbientConditions ambientConditions constrainedby TransiEnt.Components.Boundaries.Ambient.AmbientConditions "Click book icon, to edit" annotation (
    choicesAllMatching=true,
    Placement(transformation(extent={{-8,-8},{12,12}})),
    Dialog(tab="Ambience", group="Varying ambient conditions"));

  // _____________________________________________
  //
  //                   Parameters
  // _____________________________________________

  // ===== General ====

  parameter Boolean isExpertmode = true "False, show only basic parameters and variables" annotation (choices(__Dymola_checkBox=true));
  parameter Real k_H2_fraction=0.10 "Fuel fraction of hydrogen mixed with natural gas in gas turbine (Q_flow_H2/Q_flow_methane)";
  parameter Boolean isLeapYear = false "true if the observed year is a leap year" annotation (choices(__Dymola_checkBox=true));
  final parameter Modelica.Units.SI.Time lengthOfAYear=if isLeapYear then 31622400 else 31536000 "Length of one year";

  // ===== Ambience ====
  parameter Modelica.Units.SI.Pressure p_amb_const=1.013e5 "Ambient pressure" annotation (Dialog(tab="Ambience", group="Ambience parameters"));//Hamburg average 2012 (DWD, monthly average data, 11m)
  parameter Modelica.Units.SI.Temperature T_amb_const=282.48 "Ambient temperature" annotation (Dialog(tab="Ambience", group="Ambience parameters"));//Hamburg average 2012 (DWD, monthly average data, 11m)
  parameter Modelica.Units.SI.Temperature T_ground=282.48 "|Ambience|Ambience parameters|Ground temperature";
                                                                                                 //same as T_amb_const in average
  parameter Boolean variable_T_ground = false "Use variable temperature profile"
                                                                                annotation (choicesAllMatching=true,Dialog(tab="Ambience",group="Ambience parameters"));
  replaceable model Ground_Temperature =
     TransiEnt.Basics.Tables.Ambient.UndergroundTemperature_Duesseldorf_1m_3600s_TMY  constrainedby TransiEnt.Components.Boundaries.Ambient.Base.PartialTemperature
                                               "Profile for the ground temperature" annotation (choicesAllMatching=true,Dialog(tab="Ambience",group="Ambience parameters"));
  Ground_Temperature Variable_Ground_Temperature;

  parameter Real lambda=10 "degree of longitude of location" annotation(Dialog(tab="Ambience", group="Location parameters"));
  parameter Real phi=53.63 "degree of latitude of location" annotation(Dialog(tab="Ambience", group="Location parameters"));
  parameter Real timezone=1 "timezone of location (UTC+) (for Hamburg timezone=1)" annotation(Dialog(tab="Ambience", group="Location parameters"));

  // ==== Media ===


  // ==== Electric grid ====

  parameter Boolean integrateElPower=false "true if electric powers shall be integrated" annotation (Dialog(tab="Electric Grid", group="General"));
  parameter Modelica.Units.SI.Frequency f_n=50 annotation (Dialog(tab="Electric Grid", group="Nominal values (Top Level)"));
  parameter Modelica.Units.SI.Frequency delta_f_max=0.2 annotation (Dialog(tab="Electric Grid", group="Nominal values (Top Level)"));
  parameter Modelica.Units.SI.Frequency delta_f_deadband=0.01 annotation (Dialog(tab="Electric Grid", group="Nominal values (Top Level)"));
  parameter Modelica.Units.SI.Time t_SB_act=5*60 "Activation time for secondary balancing" annotation (Dialog(tab="Electric Grid", group="Nominal values (Top Level)"));
  parameter Modelica.Units.SI.Voltage v_n=110e3 annotation (Dialog(tab="Electric Grid", group="Nominal values (Top Level)"));
  parameter Modelica.Units.SI.Power P_n_low=150e9 "Nominal power of total grid at low-load" annotation (Dialog(tab="Electric Grid", group="Nominal values (Top Level)"));
  parameter Modelica.Units.SI.Power P_n_high=300e9 "Nominal power of total grid at high-load" annotation (Dialog(tab="Electric Grid", group="Nominal values (Top Level)"));
  parameter Modelica.Units.SI.Time T_grid=12 "Mechanical time constant of rotating masses in the grid" annotation (Dialog(tab="Electric Grid", group="Nominal values (Top Level)"));
  parameter Integer n_consumers=10 "Number of globaly parameterized consumers" annotation(Dialog(tab="Electric Grid", group="Optional"));
  parameter Modelica.Units.SI.ActivePower P_consumer[n_consumers]=zeros(n_consumers) "Globaly defined consumer data" annotation (Dialog(tab="Electric Grid", group="Optional"));

  replaceable TransiEnt.Grid.Electrical.Base.ExampleGenerationPark generationPark constrainedby TransiEnt.Grid.Electrical.Base.PartialGenerationPark "Properties of generaton park" annotation (Dialog(tab="Electric Grid", group="Optional"), choicesAllMatching=true);
  parameter Modelica.Units.SI.Power P_n_ref_1=generationPark.P_total "Reference power of subgrid 1, i.e. detailed grid  (for inertia constant calculation)" annotation (Dialog(tab="Electric Grid", group="Optional"));
  parameter Modelica.Units.SI.Power P_n_ref_2=P_n_high "Reference power of subgrid 2, i.e. surrounding grid  (for inertia constant calculation)" annotation (Dialog(tab="Electric Grid", group="Optional"));

  constant Integer iDetailedGrid = 1 "Index of detailed grid (convention)";
  constant Integer iSurroundingGrid = 2 "Index of detailed grid (convention)";
  parameter Modelica.Units.SI.Power P_peak_1=2.2e9 "Peak load of subgrid 1, i.e. detailed grid (used for stochastic grid error models)" annotation (Dialog(tab="Electric Grid", group="Optional"));
  parameter Modelica.Units.SI.Power P_peak_2=P_n_high "Peak load of subgrid 2, i.e. surrounding grid (used for stochastic grid error models)" annotation (Dialog(tab="Electric Grid", group="Optional"));
  //parameter Modelica.SIunits.Time T_A=5 "Time constant of all rotating masses" annotation (Dialog(tab="Electric Grid", group="Nominal values (Top Level)"));

  parameter Boolean use_reference_polar_wheel = false "If true, polar wheel angle is used as voltage angle refernce. Set true when using more than one frequency." annotation (Dialog(tab="Electric Grid", group="Optional"));

  parameter Boolean idealSuperstructLocalGrid = false "Enable to centralize ideal voltage control in each superstructure" annotation(Dialog(tab="Electric Grid", group="Optional"));
  parameter Real cosphi=1 "Reactive power factor"  annotation(Dialog(tab="Electric Grid", group="Optional"));

  //Modelica.SIunits.Frequency f_global(start=50) "global Frequency in the electric grid" annotation (Dialog(tab="Electric Grid", group="Nominal values (Top Level)"));


  // ==== Gas grid ====

  replaceable parameter TransiEnt.Basics.Media.Gases.VLE_VDIWA_NG7_H2_SRK_var gasModel1 constrainedby TILMedia.VLEFluidTypes.BaseVLEFluid "Medium name of real gas model" annotation (Dialog(tab="Media and Materials", group="TransiEnt-based models: Gas Grid"), choicesAllMatching);
  replaceable parameter TransiEnt.Basics.Media.Gases.Gas_VDIWA_NG7_H2_var gasModel2 constrainedby TILMedia.GasTypes.BaseGas "Medium name of ideal gas model" annotation (Dialog(tab="Media and Materials", group="TransiEnt-based models: Gas Grid"), choicesAllMatching);
  replaceable parameter TransiEnt.Basics.Media.Gases.VLE_VDIWA_H2_SRK gasModel3 constrainedby TILMedia.VLEFluidTypes.BaseVLEFluid "Medium name of real gas model" annotation (Dialog(tab="Media and Materials", group="TransiEnt-based models: Gas Grid"), choicesAllMatching);
  replaceable parameter TransiEnt.Basics.Media.Gases.VLE_VDIWA_SG4_var gasModel4 constrainedby TILMedia.VLEFluidTypes.BaseVLEFluid "Medium name of real gas model" annotation (Dialog(tab="Media and Materials", group="TransiEnt-based models: Gas Grid"), choicesAllMatching);
  replaceable parameter TransiEnt.Basics.Media.Gases.Gas_ExhaustGas exhaustGasModel constrainedby TILMedia.GasTypes.BaseGas "Medium name of ideal exhaust gas model" annotation (Dialog(tab="Media and Materials", group="TransiEnt-based models: Gas Grid"), choicesAllMatching);

  //replaceable parameter TransiEnt.Producer.Combined.CHPPackage.Records.ExhaustGas_VLEType exhaustVLEModel
  // constrainedby TILMedia.VLEFluidTypes.BaseVLEFluid "Medium name of real exhaust gas model" annotation(Dialog(tab="Media and Materials", group="TransiEnt-based models: Gas Grid"),choicesAllMatching);

  parameter Modelica.Units.SI.Pressure p_eff_1=25e2 "|Gas Grid|Nominal Values|Effective gauge pressure at distribution level";
  parameter Modelica.Units.SI.Pressure p_eff_2=16e5 "|Gas Grid|Nominal Values|Effective gauge pressure at distribution level";
  parameter Modelica.Units.SI.Pressure p_eff_3=25e5 "|Gas Grid|Nominal Values|Effective gauge pressure at distribution level";

  parameter Modelica.Units.SI.VolumeFraction phi_H2max=0.1 "|Gas Grid|Parameters|Maximum admissible volume fraction of H2 in NGH2 at STP";
  parameter Real f_gasDemand=2.7097280217 "|Gas Grid|Parameters|Scaling factor gas demand"; //3.0104424893;//

  parameter Boolean useConstCompInGasComp=false "true if equations for constant gas composition should be used in gas components" annotation (Dialog(tab="Gas Grid", group="Parameters"));
  parameter Integer initOptionGasPipes=0 "Type of initialization" annotation (Dialog(tab="Gas Grid", group="Parameters"), choices(
      choice=0 "Use guess values",
      choice=208 "Steady pressure and enthalpy",
      choice=201 "Steady pressure",
      choice=202 "Steady enthalpy",
      choice=210 "Steady density"));
  parameter Integer massBalanceGasPipes=1 "Mass balance and species balance fomulation" annotation (Dialog(tab="Gas Grid", group="Parameters"), choices(
      choice=1 "ClaRa formulation",
      choice=2 "TransiEnt formulation 1a",
      choice=3 "TransiEnt formulation 1b"));
  parameter Integer variableCompositionEntriesGasPipes[:](min=0)={0} "Entries of medium vector in gas pipes which are supposed to be completely variable" annotation(Dialog(tab="Gas Grid",group="Parameters",enable=not useConstCompInGasComp));

  parameter Modelica.Units.SI.Height roughnessGasPipes=0.1e-3 "Absolute roughness of gas pipes" annotation (Dialog(tab="Gas Grid", group="Parameters"));

  // ===== Expert Settings ====

  parameter Real Td = 1e-3 "Time constant of derivative calculation" annotation (Dialog(tab="Expert Settings", group="Numerical Properties"));
  parameter Boolean useThresh = false "Use threshold in gradient limiters" annotation (Dialog(tab="Expert Settings", group="Numerical Properties"));
  parameter Real thres = 1e-7 "If abs(u-y)< thres, y becomes a simple pass through of u. Increasing thres can improve simulation speed. However to large values can make the simulation unstable. 
     A good starting point is the choice thres = tolerance/1000." annotation (Dialog(tab="Expert Settings", group="Numerical Properties"));

  parameter Modelica.Units.SI.MassFlowRate m_flow_small(min=0) = 1e-2 "Default small mass flow rate for regularization of laminar and zero flow" annotation (Dialog(tab="Expert Settings", group="Numerical Properties"));

  parameter Modelica.Units.SI.Pressure p_small(min=0) = 1e3 "Default small pressure e.g. used for error handling in compressors" annotation (Dialog(tab="Expert Settings", group="Numerical Properties"));

  parameter Modelica.Units.SI.Velocity v_wind_small(min=0) = 0.1 "Default small wind velocity for startup of wind turbines (neglectable wind power)" annotation (Dialog(tab="Expert Settings", group="Numerical Properties"));

  // ==== Table Interpolation ====

  parameter Modelica.Blocks.Types.Smoothness tableInterpolationSmoothness=Modelica.Blocks.Types.Smoothness.LinearSegments "Smoothness of table interpolation"
      annotation(Dialog(tab="Expert Settings", group="table data interpretation"));

  parameter Modelica.Units.SI.Power P_el_small=1 "Small power flow considered as zero" annotation (Dialog(tab="Expert Settings", group="Singularities"));
  parameter Modelica.Units.SI.Energy E_small=3600 "Small energy quantity, considered as zero" annotation (Dialog(tab="Expert Settings", group="Singularities"));
  parameter Modelica.Units.SI.HeatFlowRate Q_flow_small=1 "Small power flow considered as zero" annotation (Dialog(tab="Expert Settings", group="Singularities"));

  // ==== Economics and Emissions ====

  //Costs Cavern
  parameter TransiEnt.Basics.Units.Time_year lifeTime_Cavern=100 annotation (Dialog(tab="Costs", group="Plant life for annuitiy"));
  parameter TransiEnt.Basics.Units.Time_year lifeTime_Electrolyzer=20 annotation (Dialog(tab="Costs", group="Plant life for annuitiy"));
  //Costs Oil
  parameter TransiEnt.Basics.Units.MonetaryUnitPerEnergy Cfue_Oil(displayUnit="EUR/J") = 45*1/159*1/36e6 annotation (Dialog(tab="Costs", group="Fuel costs in EUR/J_fuel"));
                                                                                                                                                                            //X [EUR/barrel]*1/159[barrel/liter]*1/36[liter/MJ]*3.6e9[J/MWh]


  // //Subsidies, Feed in Tariffs and selling prices


  //Selling prices in EUR/J
  parameter TransiEnt.Basics.Units.MonetaryUnitPerEnergy sellPriceElEnergy(displayUnit="EUR/J") = 210/3.6e9 "EUR/J_el" annotation (Dialog(tab="PricesAndSubsidies", group="Sale prices"));
  parameter TransiEnt.Basics.Units.MonetaryUnitPerEnergy sellPriceDHNHeat(displayUnit="EUR/J") = 81/3.6e9 "EUR/J_th" annotation (Dialog(tab="PricesAndSubsidies", group="Sale prices"));
                                                                                                                                                                                       //Source: Bundeskartellamt, 2012: Sektoruntersuchung Fernwärme
  parameter TransiEnt.Basics.Units.MonetaryUnitPerEnergy PhelixBaseYearFuture(displayUnit="EUR/J") = 29/3.6e9 "EUR/J_el" annotation (Dialog(tab="PricesAndSubsidies", group="Sale prices"));

                                                                                                                                                                                           //Source: https://www.eex.com/en/market-data/power/futures/phelix-futures oder http://www.finanzen.net/rohstoffe/eex-strom-phelix-baseload-year-future
  //Electricity selling prices time series in EUR/kWh
  replaceable TransiEnt.Basics.Tables.ElectricGrid.ElectricityPrices.SpotPriceElectricity_Phelix_DayAhead_3600s_2011 electricityPrice constrainedby TransiEnt.Basics.Tables.ElectricGrid.ElectricityPrices.GenericElectricityPriceTable "Electricity market prices in EUR per kWh" annotation (Dialog(tab="PricesAndSubsidies"), choicesAllMatching);

  //Demand-related cost
  //Free Energy
  parameter TransiEnt.Basics.Units.MonetaryUnitPerEnergy Cspec_demAndRev_free=0 "Free energy" annotation (Dialog(tab="PricesAndSubsidies", group="Demand-related costs"));
  //Electricity
  parameter TransiEnt.Basics.Units.MonetaryUnitPerEnergy Cspec_demAndRev_el_70_150_GWh=103.3/3.6e9 "Second half of 2014 for electricity usage of 70-150 GWh/a from Eurostat http://ec.europa.eu/eurostat/web/energy/data/database" annotation (Dialog(tab="PricesAndSubsidies", group="Demand-related costs"));
  //Heat

  //Gas and Fuel
  parameter TransiEnt.Basics.Units.MonetaryUnitPerEnergy Cspec_demAndRev_gas_diesel=0.8365/(0.840*43e6) "(=83.37 EUR/MWh) 0.8365 EUR/l_diesel in 2016 (Statistisches Bundesamt https://www.destatis.de/DE/Publikationen/Thematisch/Preise/Energiepreise/EnergiepreisentwicklungPDF_5619001.pdf?__blob=publicationFile), calorific value 43 MJ/kg, density 0.840 kg/l (Forschungsstelle fuer Energiewirtschaft e.V., 2010. Basisdaten zur Bereitstellung elektrischer Energie)" annotation (Dialog(tab="PricesAndSubsidies", group="Demand-related costs"));
  //Other
  parameter Real Cspec_demAndRev_other_free=0 "Free other resource" annotation (Dialog(tab="PricesAndSubsidies", group="Demand-related costs"));
  parameter Real Cspec_demAndRev_other_water=0.00185 "EUR/m3 water, 1000kg/m3, including taxes, without waste water cost, https://www.hamburgwasser.de/privatkunden/service/gebuehren-abgaben-preise/" annotation (Dialog(tab="PricesAndSubsidies", group="Demand-related costs"));

  // _____________________________________________
  //
  //             Interfaces
  // _____________________________________________

  TransiEnt.Basics.Interfaces.General.TemperatureCelsiusOut T_amb_var=ambientConditions.temperature.value "Temperature in degC (from component ambientConditions)";
  TransiEnt.Basics.Interfaces.Ambient.VelocityOut v_wind=ambientConditions.wind.value "Wind speed (from component ambientConditions)";
  TransiEnt.Basics.Interfaces.Ambient.IrradianceOut i_global=ambientConditions.globalSolarRadiation.value "Global solar radiation (from component ambientConditions)";
  TransiEnt.Basics.Interfaces.Ambient.IrradianceOut i_direct=ambientConditions.directSolarRadiation.value "Direct solar radiation (from component ambientConditions)";
  TransiEnt.Basics.Interfaces.Ambient.IrradianceOut i_diffuse=ambientConditions.diffuseSolarRadiation.value "Diffuse solar radiation (from component ambientConditions)";
  Modelica.Blocks.Interfaces.RealOutput T_ground_var(value=if variable_T_ground then Variable_Ground_Temperature.value else T_ground)  "Diffuse solar radiation (from component ambientConditions)";

   annotation ( defaultComponentName="simCenter",
    defaultComponentPrefixes="inner",
    missingInnerMessage=
        "Your model is using an outer \"simCenter\" but it does not contain an inner \"simCenter\" component. Drag model TransiEnt.SimCenter into your model to make it work.", Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics),
    Documentation(info="<html>
<h4><span style=\"color:#008000\">1. Purpose of model</span></h4>
<p>Global parameters for all models depending TransiEnt core library and Clara library.</p>
<h4><span style=\"color:#008000\">2. Level of detail, physical effects considered, and physical insight</span></h4>
<p>(Purely technical component without physical modeling.)</p>
<h4><span style=\"color:#008000\">3. Limits of validity </span></h4>
<p>(Purely technical component without physical modeling.)</p>
<h4><span style=\"color:#008000\">4. Interfaces</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color:#008000\">5. Nomenclature</span></h4>
<p>(no elements)</p>
<h4><span style=\"color:#008000\">6. Governing Equations</span></h4>
<p>(no equations)</p>
<h4><span style=\"color:#008000\">7. Remarks for Usage</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color:#008000\">8. Validation</span></h4>
<p>(no validation or testing necessary)</p>
<h4><span style=\"color:#008000\">9. References</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color:#008000\">10. Version History</span></h4>
<p>Model created by Lisa Andresen (andresen@tuhh.de) on Mon Aug 18 2014</p>
</html>"));
end SimCenter;
