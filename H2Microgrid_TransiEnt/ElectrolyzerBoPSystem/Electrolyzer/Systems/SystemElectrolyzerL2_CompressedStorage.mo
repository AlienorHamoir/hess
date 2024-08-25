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
    redeclare model electrolyzerVoltage = TransiEnt.Producer.Gas.Electrolyzer.Base.Physics.Voltage.V_cell1,
    redeclare model electrolyzerTemperature = H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.PhysicsSubmodels.Temperature_mod,
    redeclare model electrolyzerPressures = TransiEnt.Producer.Gas.Electrolyzer.Base.Physics.Pressures.Pressures1,
    redeclare model electrolyzerMassFlow = TransiEnt.Producer.Gas.Electrolyzer.Base.Physics.MassFlow.MassFlow0thOrderDynamics) annotation (Placement(transformation(extent={{-46,-18},{-6,18}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn CompPower "Electrical power from the storage compressor" annotation (Placement(transformation(extent={{-100,-76},{-60,-38}}), iconTransformation(extent={{-100,-76},{-60,-38}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_EL_tot "Total electrolyzer system consumption " annotation (Placement(transformation(extent={{94,-42},{116,-20}}), iconTransformation(extent={{94,-42},{116,-20}})));
  TransiEnt.Producer.Gas.Electrolyzer.Controller.MinMaxController minMaxController(
    P_el_n=P_el_n,
    P_el_max=P_el_max,
    P_el_min=P_el_min)             annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={-36,56})));
  CoolingSystem.HeatPortCooling.CoolingModel_Pump coolingModel_Pump annotation (Placement(transformation(extent={{46,18},{84,56}})));
  Modelica.Blocks.Interfaces.RealInput T_env "Prescribed boundary temperature from weather file" annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={38,108}), iconTransformation(
        extent={{-14.6963,37.304},{11.3045,11.3041}},
        rotation=-90,
        origin={58.696,101.304})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_EL_aux "electrolyzer auxiliaries power consumption" annotation (Placement(transformation(extent={{62,-90},{100,-52}}), iconTransformation(extent={{94,-82},{118,-58}})));

equation

  connect(epp, electrolyzer.epp) annotation (Line(
      points={{-100,0},{-46,0}},
      color={0,135,135},
      thickness=0.5));
  connect(electrolyzer.CompressorPowerIn, CompPower) annotation (Line(points={{-50,-13.68},{-60,-13.68},{-60,-58},{-80,-58},{-80,-57}}, color={0,127,127}));
  connect(minMaxController.P_el_ely, electrolyzer.P_el_set) annotation (Line(
      points={{-36,45.2},{-35.2,45.2},{-35.2,16.92}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(minMaxController.P_el_set, P_el_set) annotation (Line(points={{-36,67},{-36,82},{0,82},{0,108}}, color={0,127,127}));
  connect(electrolyzer.electrolyzerPowerOut, P_EL_tot) annotation (Line(
      points={{-40.4,-12.6},{-40.4,-32},{54,-32},{54,-31},{105,-31}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(electrolyzer.temperatureOut, coolingModel_Pump.T_op) annotation (Line(points={{-40.4,-5.04},{-62,-5.04},{-62,24},{32,24},{32,50.3},{46,50.3}},
                                                                                                                                                     color={0,0,127}));
  connect(electrolyzer.heat, coolingModel_Pump.heatPortCooling) annotation (Line(points={{-6,-11.88},{36,-11.88},{36,23.51},{46.19,23.51}},
                                                                                                                                color={191,0,0}));
  connect(coolingModel_Pump.P_CoolPump, electrolyzer.CoolingPumpPowerIn) annotation (Line(points={{83.81,44.41},{92,44.41},{92,14},{-4,14},{-4,2.88},{-3.2,2.88}},color={0,0,127}));
  connect(coolingModel_Pump.T_env, T_env) annotation (Line(points={{46,37},{42,37},{42,38},{38,38},{38,108}}, color={0,0,127}));
  connect(electrolyzer.electrolyzerPowerOutAux, P_EL_aux) annotation (Line(
      points={{-40.4,-16.2},{-40.4,-70},{82,-70},{82,-71},{81,-71}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(gasPortOut, electrolyzer.gasPortOut) annotation (Line(
      points={{5,-94},{6,-94},{6,0},{-6,0}},
      color={255,255,0},
      thickness=1.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid), Text(
          extent={{-76,38},{76,-42}},
          textColor={255,255,255},
          textStyle={TextStyle.Bold},
          textString="Electrolyzer
L2")}),                                     Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>Electrolyzer system with L2 electrolyzer model and its accounted BoP, consisting of a cooling model with a pump. All models are found in H2Microgrid_TransiEnt and are based on existing TransiEnt models. </p>
<p>The high-pressure storage compressor power is used as input in the electrolyzer model.</p>
<p>The electrolyzer system also contains a power controller found in TransiEnt library.</p>
<p><br>Tested in the check models &quot;H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Tests.Test_PEMElectrolyzerL2_CompressedStorage&quot;</p>
<p>Simplified feed-in station model, without connection to a natural gas grid</p>
</html>"));
end SystemElectrolyzerL2_CompressedStorage;
