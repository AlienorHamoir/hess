within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.Systems;
model SystemElectrolyzerL2_CompressedStorage

  extends TransiEnt.Producer.Gas.Electrolyzer.Base.PartialFeedInStation(
    gasPortOut(Medium=medium),
    break m_flow_feedIn);


  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_coolant=simCenter.fluid1;

  parameter Modelica.Units.SI.ActivePower P_el_n=5e3 "nominal power of electrolyser" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.ActivePower P_el_max=1.0*P_el_n "Maximum power of electrolyzer" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.ActivePower P_el_min=0.1*P_el_n "Minimal power of electrolyzer" annotation (Dialog(tab="General", group="Electrolyzer"));
//  parameter Real eta_inv_n=0.956 "nominal efficiency of the inverter" annotation (Dialog(group="Fundamental Definitions")); // we neglect all losses relative to grid constraints, including inverter efficiency


  parameter Modelica.Units.SI.MassFlowRate m_flow_start=0.0 "Sets initial value for m_flow from a buffer" annotation (Dialog(tab="General", group="Initialization"));
// parameter SI.Temperature T_Init=283.15 "Sets initial value for T" annotation (Dialog(tab="General", group="Initialization"));

//// If electrolyzer can work in overload (not our case)
//   parameter Modelica.Units.SI.ActivePower P_el_overload=1.0*P_el_n "Power at which overload region begins" annotation (Dialog(tab="General", group="Electrolyzer"));
//   parameter Modelica.Units.SI.ActivePower P_el_cooldown=0.5*P_el_n "Power below which cooldown of electrolyzer starts" annotation (Dialog(group="Electrolyzer"));
//   parameter Real t_overload=0.5*3600 "Maximum admissible time the electrolyzer can work in overload region in seconds" annotation (Dialog(tab="General", group="Electrolyzer"));
//   parameter Real coolingToHeatingRatio=1 "Defines how much faster electrolyzer cools down than heats up" annotation (Dialog(tab="General", group="Electrolyzer"));
//   parameter Integer startState=1 "Initial state of the electrolyzer (1: ready to overheat, 2: working in overload, 3: cooling down)" annotation (Dialog(tab="General", group="Electrolyzer"));

  parameter Modelica.Units.SI.AbsolutePressure p_out=30e5 "Hydrogen output pressure from electrolyser" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.Temperature T_out=273.15 + 40 "Hydrogen output temperature from electrolyser" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Real specificWaterConsumption=11 "Mass of water per mass of hydrogen" annotation (Dialog(tab="General", group="Electrolyzer"));

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
        origin={105,65}), iconTransformation(
        extent={{-12,-12},{12,12}},
        rotation=0,
        origin={106,-16})));
  PEMElectrolyzer_L2 electrolyzer(
    useFluidCoolantPort=false,
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
    medium=medium,
    eta_GCV_EL(start=0),
    redeclare model electrolyzerVoltage = H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.PhysicsSubmodels.V_cell,
    redeclare model electrolyzerTemperature = H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.PhysicsSubmodels.Temperature_mod,
    redeclare model electrolyzerPressures = TransiEnt.Producer.Gas.Electrolyzer.Base.Physics.Pressures.Pressures1,
    redeclare model electrolyzerMassFlow = TransiEnt.Producer.Gas.Electrolyzer.Base.Physics.MassFlow.MassFlow0thOrderDynamics) annotation (Placement(transformation(extent={{-66,-18},{-26,18}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn CompressorPower "Electrical power from the storage compressor" annotation (Placement(transformation(extent={{-120,-72},{-90,-44}}), iconTransformation(extent={{-120,-72},{-90,-44}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_electrolyzer_tot annotation (Placement(transformation(extent={{94,-82},{118,-58}}), iconTransformation(extent={{94,-82},{118,-58}})));
  TransiEnt.Producer.Gas.Electrolyzer.Controller.MinMaxController minMaxController(
    P_el_n=P_el_n,
    P_el_max=P_el_max,
    P_el_min=P_el_min)             annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={-54,46})));
  CoolingSubystem.HeatPortCooling.CoolingModel_Valve coolingModel_Valve annotation (Placement(transformation(extent={{44,-38},{78,-12}})));
  Modelica.Blocks.Sources.Constant nulPower(k=0) annotation (Placement(transformation(extent={{-24,22},{-10,36}})));
protected
  TransiEnt.Components.Sensors.RealGas.MassFlowSensor massflowSensor_ely(medium=medium, xiNumber=0)         annotation (Placement(transformation(
        extent={{7,6},{-7,-6}},
        rotation=180,
        origin={-6,6})));

equation
  connect(massflowSensor_ely.m_flow, H2massFlowRateOutElectrolyzer) annotation (Line(points={{1.7,6},{88,6},{88,65},{105,65}},
                                                                                                                        color={0,0,127}));

  connect(gasPortOut, massflowSensor_ely.gasPortOut) annotation (Line(
      points={{5,-94},{1,-94},{1,0}},
      color={255,255,0},
      thickness=1.5));
  connect(epp, electrolyzer.epp) annotation (Line(
      points={{-100,0},{-66,0}},
      color={0,135,135},
      thickness=0.5));
  connect(electrolyzer.gasPortOut, massflowSensor_ely.gasPortIn) annotation (Line(
      points={{-26,0},{-13,0}},
      color={255,255,0},
      thickness=1.5));
  connect(electrolyzer.storageCompressorPowerIn, CompressorPower) annotation (Line(points={{-70,-13.68},{-84,-13.68},{-84,-58},{-105,-58}}, color={0,127,127}));
  connect(minMaxController.P_el_ely, electrolyzer.P_el_set) annotation (Line(
      points={{-54,35.2},{-53.6,36},{-53.6,21.6}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(minMaxController.P_el_set, P_el_set) annotation (Line(points={{-54,57},{-54,82},{0,82},{0,108}}, color={0,127,127}));
  connect(electrolyzer.electrolyzerPowerOut, P_electrolyzer_tot) annotation (Line(
      points={{-57.2,-12.6},{-57.2,-70},{106,-70}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(coolingModel_Valve.T_op, electrolyzer.temperatureOut) annotation (Line(points={{41.96,-16.42},{-86,-16.42},{-86,-5.76},{-57.2,-5.76}}, color={0,0,127}));
  connect(electrolyzer.heat, coolingModel_Valve.heatPortCooling) annotation (Line(points={{-26,-11.88},{34,-11.88},{34,-42},{44,-42},{44,-36.96}}, color={191,0,0}));
  connect(nulPower.y, electrolyzer.coolingPumpPowerIn) annotation (Line(points={{-9.3,29},{-6,29},{-6,16},{-20,16},{-20,12},{-22,12},{-22,3.6}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid), Text(
          extent={{-76,38},{76,-42}},
          textColor={255,255,255},
          textStyle={TextStyle.Bold},
          textString="Electrolyzer
L2")}),                                     Diagram(coordinateSystem(preserveAspectRatio=false)));
end SystemElectrolyzerL2_CompressedStorage;
