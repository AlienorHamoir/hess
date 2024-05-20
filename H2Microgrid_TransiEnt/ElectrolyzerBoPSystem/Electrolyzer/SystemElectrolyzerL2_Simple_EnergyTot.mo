within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer;
model SystemElectrolyzerL2_Simple_EnergyTot

  extends TransiEnt.Producer.Gas.Electrolyzer.Base.PartialFeedInStation(gasPortOut(Medium=medium));

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;


  parameter Modelica.Units.SI.ActivePower P_el_n=5.5e3 "nominal power of electrolyser" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.ActivePower P_el_max=2.0*P_el_n "Maximum power of electrolyzer" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.ActivePower P_el_min=0.05*P_el_n "Minimal power of electrolyzer" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.ActivePower P_el_overload=1.0*P_el_n "Power at which overload region begins" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.ActivePower P_el_cooldown=0.5*P_el_n "Power below which cooldown of electrolyzer starts" annotation (Dialog(group="Electrolyzer"));
  parameter Modelica.Units.SI.MassFlowRate m_flow_start=0.0 "Sets initial value for m_flow from a buffer" annotation (Dialog(tab="General", group="Initialization"));
// parameter SI.Temperature T_Init=283.15 "Sets initial value for T" annotation (Dialog(tab="General", group="Initialization"));
// parameter Modelica.Units.SI.Efficiency eta_n(
//    min=0,
//    max=1)=0.75 "Nominal efficency refering to the GCV (min = 0, max = 1)" annotation (Dialog(tab="General", group="Electrolyzer"));
//  parameter Modelica.Units.SI.Efficiency eta_scale(
//    min=0,
//    max=1)=0 "Sets a with increasing input power linear degrading efficiency coefficient (min = 0, max = 1)" annotation (Dialog(tab="General", group="Electrolyzer")); for electrolyzer L1 but not L2

  parameter Real t_overload=0.5*3600 "Maximum admissible time the electrolyzer can work in overload region in seconds" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Real coolingToHeatingRatio=1 "Defines how much faster electrolyzer cools down than heats up" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Integer startState=1 "Initial state of the electrolyzer (1: ready to overheat, 2: working in overload, 3: cooling down)" annotation (Dialog(tab="General", group="Electrolyzer"));

  parameter Modelica.Units.SI.AbsolutePressure p_out=30e5 "Hydrogen output pressure from electrolyser" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.Temperature T_out=273.15 + 10 "Hydrogen output temperature from electrolyser" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Real specificWaterConsumption=10 "Mass of water per mass of hydrogen" annotation (Dialog(tab="General", group="Electrolyzer"));

  parameter Real k=1e8 "Gain for feed-in control" annotation (Dialog(tab="General", group="Controller"));
  parameter Real Ti=0.1 "Ti for feed-in control" annotation (Dialog(tab="General", group="Controller"));
  parameter Real Td=0.1 "Td for feed-in control" annotation (Dialog(tab="General", group="Controller"));

  parameter Boolean useFluidCoolantPort=false "choose if fluid port for coolant shall be used" annotation (Dialog(enable=not useHeatPort,group="Coolant"));
  parameter Boolean useHeatPort=false "choose if heat port for coolant shall be used" annotation (Dialog(enable=not useFluidCoolantPort,group="Coolant"));
  parameter Boolean useVariableCoolantOutputTemperature=true "choose if temperature of cooland output shall be defined by input" annotation (Dialog(enable=useFluidCoolantPort,group="Coolant"));
  parameter Modelica.Units.SI.Temperature T_out_coolant_target=50 + 273.15 "output temperature of coolant - will be limited by temperature which is technically feasible" annotation (Dialog(enable=useFluidCoolantPort, group="Coolant"));
  parameter Boolean externalMassFlowControl=false "choose if heat port for coolant shall be used" annotation (Dialog(enable=useFluidCoolantPort and (not useVariableCoolantOutputTemperature), group="Coolant"));

//   replaceable model Charline = TransiEnt.Producer.Gas.Electrolyzer.Base.ElectrolyzerEfficiencyCharlineSilyzer200        constrainedby TransiEnt.Producer.Gas.Electrolyzer.Base.PartialElectrolyzerEfficiencyCharline        "Calculate the efficiency" annotation (Dialog(group="Electrolyzer"),__Dymola_choicesAllMatching=true);
//   replaceable model Dynamics = TransiEnt.Producer.Gas.Electrolyzer.Base.ElectrolyzerDynamics0thOrder        constrainedby TransiEnt.Producer.Gas.Electrolyzer.Base.PartialElectrolyzerDynamics        "Dynamic behavior of electrolyser" annotation (Dialog(group="Electrolyzer"),__Dymola_choicesAllMatching=true);
//   replaceable model PressureLossAtOutlet = ClaRa.Components.VolumesValvesFittings.Fittings.Fundamentals.NoFriction constrainedby ClaRa.Components.VolumesValvesFittings.Fittings.Fundamentals.BaseDp "Pressure loss at outlet (can help reducing the system of equations)" annotation(choicesAllMatching=true);

model ElectrolyzerRecord
    extends TransiEnt.Basics.Icons.Record;
    input SI.Power P_el "Consumed electric power";
    input SI.Energy W_el "Consumed electric energy";
    input SI.Mass mass_H2 "Produced hydrogen mass";
    input SI.Efficiency eta_NCV "Electroyzer efficiency based on NCV";
    input SI.Efficiency eta_GCV "Electroyzer efficiency based on GCV";
end ElectrolyzerRecord;

//     final eta_n=eta_n,
//     final eta_scale=eta_scale,

//     redeclare model Dynamics = Dynamics,
//     redeclare model Charline = Charline,
//     redeclare model CostSpecsGeneral = CostSpecsElectrolyzer,
  TransiEnt.Basics.Interfaces.General.MassFlowRateOut H2massFlowRateOutElectrolyzer annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={50,110}), iconTransformation(
        extent={{-12,-12},{12,12}},
        rotation=0,
        origin={106,-6})));
  TransiEnt.Basics.Interfaces.General.TemperatureIn T_set_coolant_out if useVariableCoolantOutputTemperature annotation (Placement(transformation(extent={{112,26},{88,50}}), iconTransformation(extent={{112,26},{88,50}})));
  TransiEnt.Basics.Interfaces.Thermal.FluidPortOut fluidPortOut(Medium=simCenter.fluid1) if useFluidCoolantPort annotation (Placement(transformation(extent={{-70,-108},{-50,-88}}), iconTransformation(extent={{-70,-108},{-50,-88}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heat if useHeatPort annotation (Placement(transformation(extent={{90,-80},{114,-56}}), iconTransformation(extent={{90,-80},{114,-56}})));
  TransiEnt.Basics.Interfaces.Thermal.FluidPortIn fluidPortIn(Medium=simCenter.fluid1) if useFluidCoolantPort annotation (Placement(transformation(extent={{50,-108},{70,-88}}), iconTransformation(extent={{50,-108},{70,-88}})));
  PEMElectrolyzer_L2_EnergyTot electrolyzer1(
    useFluidCoolantPort=useFluidCoolantPort,
    useHeatPort=true,
    externalMassFlowControl=false,
    useVariableCoolantOutputTemperature=false,
    T_out_coolant_target=T_out_coolant_target,
    usePowerPort=true,
    P_el_n=P_el_n,
    P_el_max=P_el_max,
    integrateH2Flow=true,
    integrateElPower=true,
    T_out=T_out,
    medium=medium) annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
  TransiEnt.Basics.Interfaces.General.TemperatureOut temperatureOutElectrolyzer annotation (Placement(transformation(extent={{94,-2},{118,22}}), iconTransformation(extent={{94,-2},{118,22}})));
protected
  TransiEnt.Components.Sensors.RealGas.MassFlowSensor massflowSensor_ely(medium=medium, xiNumber=0)         annotation (Placement(transformation(
        extent={{7,6},{-7,-6}},
        rotation=180,
        origin={-6,6})));

equation
  connect(massflowSensor_ely.m_flow, H2massFlowRateOutElectrolyzer) annotation (Line(points={{1.7,6},{6,6},{6,84},{50,84},{50,110}},
                                                                                                                        color={0,0,127}));

//       if integrateMassFlow then
//   der(mass_H2_fedIn)=-gasPortOut.m_flow;
// else
//   mass_H2_fedIn=0;
// end if;

  if usePowerPort then
  end if;

  if useFluidCoolantPort then
  end if;
  if useHeatPort then
    connect(electrolyzer1.heat, heat) annotation (Line(points={{-40,-6.6},{84,-6.6},{84,-68},{102,-68}}, color={191,0,0}));
  end if;
  if useVariableCoolantOutputTemperature then
  end if;

  connect(gasPortOut, massflowSensor_ely.gasPortOut) annotation (Line(
      points={{5,-94},{1,-94},{1,0}},
      color={255,255,0},
      thickness=1.5));
  connect(electrolyzer1.epp, epp) annotation (Line(
      points={{-60,0},{-100,0}},
      color={0,135,135},
      thickness=0.5));
  connect(electrolyzer1.gasPortOut, massflowSensor_ely.gasPortIn) annotation (Line(
      points={{-40,0},{-13,0}},
      color={255,255,0},
      thickness=1.5));
  connect(P_el_set, electrolyzer1.P_el_set) annotation (Line(points={{0,108},{0,26},{-53.8,26},{-53.8,12}},                      color={0,127,127}));
  connect(electrolyzer1.heat, heat) annotation (Line(points={{-40,-6.6},{84,-6.6},{84,-68},{102,-68}}, color={191,0,0}));
  connect(temperatureOutElectrolyzer, electrolyzer1.temperatureOut) annotation (Line(points={{106,10},{8,10},{8,-8},{-36,-8},{-36,-14},{-64,-14},{-64,-3.2},{-55.6,-3.2}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio=false)));
end SystemElectrolyzerL2_Simple_EnergyTot;
