within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer;
model SystemElectrolyzerL2_Simple_EnergyTot_Inv_Cooling

  extends TransiEnt.Producer.Gas.Electrolyzer.Base.PartialFeedInStation(gasPortOut(Medium=medium), break m_flow_feedIn);

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_coolant=simCenter.fluid1;

  parameter Modelica.Units.SI.ActivePower P_el_n=5.5e3 "nominal power of electrolyser" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.ActivePower P_el_max=2.0*P_el_n "Maximum power of electrolyzer" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.ActivePower P_el_min=0.05*P_el_n "Minimal power of electrolyzer" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.ActivePower P_el_overload=1.0*P_el_n "Power at which overload region begins" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.ActivePower P_el_cooldown=0.5*P_el_n "Power below which cooldown of electrolyzer starts" annotation (Dialog(group="Electrolyzer"));
  parameter Modelica.Units.SI.MassFlowRate m_flow_start=0.0 "Sets initial value for m_flow from a buffer" annotation (Dialog(tab="General", group="Initialization"));
// parameter SI.Temperature T_Init=283.15 "Sets initial value for T" annotation (Dialog(tab="General", group="Initialization"));

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
  parameter Boolean useHeatPort=true "choose if heat port for coolant shall be used" annotation (Dialog(enable=not useFluidCoolantPort,group="Coolant"));
  parameter Boolean useVariableCoolantOutputTemperature=false "choose if temperature of cooland output shall be defined by input" annotation (Dialog(enable=useFluidCoolantPort,group="Coolant"));
  parameter Modelica.Units.SI.Temperature T_out_coolant_target=50 + 273.15 "output temperature of coolant - will be limited by temperature which is technically feasible" annotation (Dialog(enable=useFluidCoolantPort, group="Coolant"));
  parameter Boolean externalMassFlowControl=false "choose if mass flow control for coolant shall be used" annotation (Dialog(enable=useFluidCoolantPort and (not useVariableCoolantOutputTemperature), group="Coolant"));

model ElectrolyzerRecord
    extends TransiEnt.Basics.Icons.Record;
    input SI.Power P_el "Consumed electric power";
    input SI.Energy W_el "Consumed electric energy";
    input SI.Mass mass_H2 "Produced hydrogen mass";
    input SI.Efficiency eta_NCV "Electroyzer efficiency based on NCV";
    input SI.Efficiency eta_GCV "Electroyzer efficiency based on GCV";
end ElectrolyzerRecord;

  TransiEnt.Basics.Interfaces.General.MassFlowRateOut H2massFlowRateOutElectrolyzer annotation (Placement(transformation(
        extent={{-11,-11},{11,11}},
        rotation=0,
        origin={105,7}),  iconTransformation(
        extent={{-12,-12},{12,12}},
        rotation=0,
        origin={106,-16})));
  TransiEnt.Basics.Interfaces.General.TemperatureIn T_set_coolant_out if useVariableCoolantOutputTemperature annotation (Placement(transformation(extent={{112,58},{88,82}}), iconTransformation(extent={{112,58},{88,82}})));
  TransiEnt.Basics.Interfaces.General.TemperatureOut temperatureOutElectrolyzer annotation (Placement(transformation(extent={{94,26},{118,50}}), iconTransformation(extent={{94,26},{118,50}})));
  PEMElectrolyzer_L2_EnergyTot_Inv electrolyzer(
    useFluidCoolantPort=false,
    externalMassFlowControl=false,
    useVariableCoolantOutputTemperature=false,
    T_out_coolant_target=T_out_coolant_target,
    usePowerPort=true,
    P_el_n=P_el_n,
    P_el_max=P_el_max,
    integrateH2Flow=true,
    integrateElPower=true,
    T_out=T_out,
    medium=medium,
    redeclare model electrolyzerTemperature = H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.Temperature_modPID) annotation (Placement(transformation(extent={{-66,-18},{-26,18}})));
protected
  TransiEnt.Components.Sensors.RealGas.MassFlowSensor massflowSensor_ely(medium=medium, xiNumber=0)         annotation (Placement(transformation(
        extent={{7,6},{-7,-6}},
        rotation=180,
        origin={-6,6})));

equation
  connect(massflowSensor_ely.m_flow, H2massFlowRateOutElectrolyzer) annotation (Line(points={{1.7,6},{54,6},{54,7},{105,7}},
                                                                                                                        color={0,0,127}));

  connect(gasPortOut, massflowSensor_ely.gasPortOut) annotation (Line(
      points={{5,-94},{1,-94},{1,0}},
      color={255,255,0},
      thickness=1.5));
  connect(epp, electrolyzer.epp) annotation (Line(
      points={{-100,0},{-66,0}},
      color={0,135,135},
      thickness=0.5));
  connect(electrolyzer.P_el_set, P_el_set) annotation (Line(points={{-53.6,21.6},{-53.6,82},{0,82},{0,108}}, color={0,0,127}));
  connect(electrolyzer.gasPortOut, massflowSensor_ely.gasPortIn) annotation (Line(
      points={{-26,0},{-13,0}},
      color={255,255,0},
      thickness=1.5));
  connect(electrolyzer.temperatureOut, temperatureOutElectrolyzer) annotation (Line(points={{-57.2,-5.76},{-64,-5.76},{-64,-6},{-20,-6},{-20,38},{106,38}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio=false)));
end SystemElectrolyzerL2_Simple_EnergyTot_Inv_Cooling;
