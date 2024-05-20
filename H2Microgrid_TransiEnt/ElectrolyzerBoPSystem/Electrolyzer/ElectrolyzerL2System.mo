within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer;
model ElectrolyzerL2System "Contains electrolyzer stack, H2 cooling system"

    // _____________________________________________
  //
  //          Imports and Class Hierarchy
  // _____________________________________________

  extends TransiEnt.Producer.Gas.Electrolyzer.Base.PartialFeedInStation(gasPortOut(Medium=medium));

  // _____________________________________________
  //
  //           Constants and Parameters
  // _____________________________________________

  replaceable parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel4 "Hydrogen model to be used";
  parameter Modelica.Units.SI.ActivePower P_el_n=1e6 "nominal power of electrolyser" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.ActivePower P_el_max=1.68*P_el_n "Maximum power of electrolyzer" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.ActivePower P_el_min=0.05*P_el_n "Minimal power of electrolyzer" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.ActivePower P_el_overload=1.0*P_el_n "Power at which overload region begins" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.ActivePower P_el_cooldown=P_el_n "Power below which cooldown of electrolyzer starts" annotation (Dialog(group="Electrolyzer"));
  parameter Modelica.Units.SI.MassFlowRate m_flow_start=0.0 "Sets initial value for m_flow from a buffer" annotation (Dialog(tab="General", group="Initialization"));
  //parameter SI.Temperature T_Init=283.15 "Sets initial value for T" annotation (Dialog(tab="General", group="Initialization"));
// parameter Modelica.Units.SI.Efficiency eta_n(
//    min=0,
//    max=1)=0.75 "Nominal efficency refering to the GCV (min = 0, max = 1)" annotation (Dialog(tab="General", group="Electrolyzer"));
//  parameter Modelica.Units.SI.Efficiency eta_scale(
//    min=0,
//    max=1)=0 "Sets a with increasing input power linear degrading efficiency coefficient (min = 0, max = 1)" annotation (Dialog(tab="General", group="Electrolyzer")); for electrolyzer L1 but not L2

  parameter Real k=1e8 "Gain for feed-in control" annotation (Dialog(tab="General", group="Controller"));
  parameter Real Ti=0.1 "Ti for feed-in control";
  parameter Real Td=0.1 "Td for feed-in control";

  parameter Real t_overload=0.5*3600 "Maximum admissible time the electrolyzer can work in overload region in seconds" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Real coolingToHeatingRatio=1 "Defines how much faster electrolyzer cools down than heats up" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Integer startState=1 "Initial state of the electrolyzer (1: ready to overheat, 2: working in overload, 3: cooling down)" annotation (Dialog(tab="General", group="Electrolyzer"));

  parameter Modelica.Units.SI.AbsolutePressure p_out=35e5 "Hydrogen output pressure from electrolyser" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.Temperature T_out=283.15 "Hydrogen output temperature from electrolyser" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Real specificWaterConsumption=10 "Mass of water per mass of hydrogen" annotation (Dialog(tab="General", group="Electrolyzer"));


  //replaceable model Charline = TransiEnt.Producer.Gas.Electrolyzer.Base.ElectrolyzerEfficiencyCharlineSilyzer200        constrainedby TransiEnt.Producer.Gas.Electrolyzer.Base.PartialElectrolyzerEfficiencyCharline        "Calculate the efficiency" annotation (Dialog(group="Electrolyzer"),__Dymola_choicesAllMatching=true);
  //replaceable model Dynamics = TransiEnt.Producer.Gas.Electrolyzer.Base.ElectrolyzerDynamics0thOrder        constrainedby TransiEnt.Producer.Gas.Electrolyzer.Base.PartialElectrolyzerDynamics        "Dynamic behavior of electrolyser" annotation (Dialog(group="Electrolyzer"),__Dymola_choicesAllMatching=true);
  replaceable model PressureLossAtOutlet = ClaRa.Components.VolumesValvesFittings.Fittings.Fundamentals.NoFriction constrainedby ClaRa.Components.VolumesValvesFittings.Fittings.Fundamentals.BaseDp "Pressure loss at outlet (can help reducing the system of equations)" annotation(choicesAllMatching=true);

  parameter Boolean useFluidCoolantPort=true "choose if fluid port for coolant shall be used" annotation (Dialog(enable=not useHeatPort,group="Coolant"));
  parameter Boolean useHeatPort=true "choose if heat port for coolant shall be used" annotation (Dialog(enable=not useFluidCoolantPort,group="Coolant"));
  parameter Boolean useVariableCoolantOutputTemperature=true "choose if temperature of cooland output shall be defined by input" annotation (Dialog(enable=useFluidCoolantPort,group="Coolant"));
  parameter Modelica.Units.SI.Temperature T_out_coolant_target=500 + 273.15 "output temperature of coolant - will be limited by temperature which is technically feasible" annotation (Dialog(enable=useFluidCoolantPort, group="Coolant"));
  parameter Boolean externalMassFlowControl=true "choose if heat port for coolant shall be used" annotation (Dialog(enable=useFluidCoolantPort and (not useVariableCoolantOutputTemperature), group="Coolant"));

 // _____________________________________________
  //
  //                  Variables
  // _____________________________________________
 // Modelica.Units.SI.Mass mass_H2_fedIn "Fed in hydrogen mass";
  //Modelica.Units.SI.ActivePower P_el_max_is=if controlTotalElyStorage.overloadController.state == 3 then controlTotalElyStorage.overloadController.P_el_cooldown elseif controlTotalElyStorage.feedInStorageController.storageFull then controlTotalElyStorage.feedInStorageController.limPID.y else P_el_max;

  TransiEnt.Basics.Interfaces.Thermal.FluidPortOut fluidPortOut(Medium=simCenter.fluid1) if useFluidCoolantPort annotation (Placement(transformation(extent={{90,-14},{110,6}})));
  TransiEnt.Basics.Interfaces.Thermal.FluidPortIn fluidPortIn(Medium=simCenter.fluid1) if useFluidCoolantPort annotation (Placement(transformation(extent={{92,-38},{112,-18}})));
  ElyStorageController elyStorageController annotation (Placement(transformation(extent={{-66,46},{-46,66}})));
  TransiEnt.Basics.Interfaces.General.PressureIn pressureIn "Pressure of hydrogen storage " annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={60,110})));

  TransiEnt.Basics.Interfaces.General.TemperatureIn T_set_coolant_out if useVariableCoolantOutputTemperature annotation (Placement(transformation(extent={{126,4},{86,44}})));

  inner TransiEnt.ModelStatistics modelStatistics annotation (Placement(transformation(extent={{-96,76},{-76,96}})));
  PEMElectrolyzer_L2 electrolyzer1(medium=medium)   annotation (Placement(transformation(extent={{-64,-10},{-44,10}})));
protected
  TransiEnt.Components.Sensors.RealGas.MassFlowSensor massflowSensor_ely(final medium=medium, xiNumber=0)         annotation (Placement(transformation(
        extent={{7,6},{-7,-6}},
        rotation=180,
        origin={-14,6})));

//// Got to adapt the controller to the ElyStorageController
// public
//   TransiEnt.Producer.Gas.Electrolyzer.Controller.TotalFeedInStorageController controlTotalElyStorage(
//     final k=k,
//     p_minLow=dp_Low,
//     final p_maxHigh=p_maxHigh,
//     final p_maxLow=p_maxLow,
//     final coolingToHeatingRatio=coolingToHeatingRatio,
//     P_el_n=P_el_n,
//     P_el_min=P_el_min,
//     P_el_max=P_el_max,
//     P_el_overload=P_el_overload,
//     t_overload=t_overload,
//     eta_n=eta_n,
//     eta_scale=eta_scale,
//     startState=startState,
//     p_minHigh=dp_High,
//     p_minLow_constantDemand=p_minLow_constantDemand,
//     m_flow_hydrogenDemand_constant=m_flow_hydrogenDemand_constant,
//     redeclare model Charline = Charline,
//     controllerType=Modelica.Blocks.Types.SimpleController.P,
//     StoreAllHydrogen=StoreAllHydrogen,
//     P_el_cooldown=P_el_cooldown,
//     Ti=Ti,
//     Td=Td);

equation
  connect(massflowSensor_ely.m_flow, elyStorageController.m_flow_ely) annotation (Line(points={{-6.3,6},{-2,6},{-2,51},{-46,51}},                  color={0,0,127}));
  connect(pressureIn, elyStorageController.p_storage) annotation (Line(points={{60,110},{60,57},{-46,57}}, color={0,0,127}));
  connect(m_flow_feedIn, elyStorageController.m_flow_feedIn) annotation (Line(points={{108,70},{-36,70},{-36,63},{-46,63}}, color={0,0,127}));
  connect(P_el_set, elyStorageController.P_el_set) annotation (Line(points={{0,108},{0,76},{-62,76},{-62,66}}, color={0,127,127}));
  connect(electrolyzer1.epp, epp) annotation (Line(
      points={{-64,0},{-100,0}},
      color={0,135,135},
      thickness=0.5));
  connect(electrolyzer1.P_el_set, elyStorageController.P_el_ely) annotation (Line(points={{-57.8,12},{-57.8,42},{-62,42},{-62,45}}, color={0,0,127}));
  connect(electrolyzer1.gasPortOut, massflowSensor_ely.gasPortIn) annotation (Line(
      points={{-44,0},{-21,0}},
      color={255,255,0},
      thickness=1.5));
  connect(massflowSensor_ely.gasPortOut, gasPortOut) annotation (Line(
      points={{-7,0},{-7,-4},{5,-4},{5,-94}},
      color={255,255,0},
      thickness=1.5));
  connect(fluidPortOut, electrolyzer1.fluidPortOut) annotation (Line(
      points={{100,-4},{8,-4},{8,-6},{-40,-6},{-40,-4},{-44,-4}},
      color={175,0,0},
      thickness=0.5));
  connect(fluidPortIn, electrolyzer1.fluidPortIn) annotation (Line(
      points={{102,-28},{-44,-28},{-44,-9}},
      color={175,0,0},
      thickness=0.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end ElectrolyzerL2System;
