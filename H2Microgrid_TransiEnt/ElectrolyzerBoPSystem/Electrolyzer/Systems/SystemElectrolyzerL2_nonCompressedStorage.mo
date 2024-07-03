within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.Systems;
model SystemElectrolyzerL2_nonCompressedStorage

  extends TransiEnt.Producer.Gas.Electrolyzer.Base.PartialFeedInStation(
    gasPortOut(Medium=medium) annotation(Placement(
      transformation(
        origin={400, 400},
        extent={{-10, -10}, {10, 10}},
        rotation=0))),
    break m_flow_feedIn,
    break epp);

     // Change these values to your desired position
                                          // Change the size if needed
                      // Change the rotation if needed

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
    usePowerPort=false,
    P_el_n=P_el_n,
    P_el_max=P_el_max,
    integrateH2Flow=true,
    integrateElPower=false,
    T_out=T_out,
    medium=medium,
    mass_H2(start=0.1),
    redeclare model electrolyzerVoltage = H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.PhysicsSubmodels.V_cell,
    redeclare model electrolyzerTemperature = H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.PhysicsSubmodels.Temperature_mod,
    redeclare model electrolyzerPressures = TransiEnt.Producer.Gas.Electrolyzer.Base.Physics.Pressures.Pressures1,
    redeclare model electrolyzerMassFlow = TransiEnt.Producer.Gas.Electrolyzer.Base.Physics.MassFlow.MassFlow2thOrderDynamics) annotation (Placement(transformation(extent={{-66,-18},{-26,18}})));
  Modelica.Blocks.Sources.Constant nulPower(k=0) annotation (Placement(transformation(extent={{-96,-76},{-82,-62}})));
  TransiEnt.Producer.Gas.Electrolyzer.Controller.MinMaxController minMaxController(
    P_el_n=P_el_n,
    P_el_max=P_el_max,
    P_el_min=P_el_min)             annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={-54,46})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut electrolyzerPowerOut annotation (Placement(transformation(extent={{94,-86},{114,-66}})));
  CoolingSystem.HeatPortCooling.CoolingModel_Valve coolingModel_Valve annotation (Placement(transformation(extent={{36,-38},{70,-12}})));
  Modelica.Blocks.Interfaces.RealInput T_environment "Prescribed boundary temperature from weather file" annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={80,108}), iconTransformation(
        extent={{-14.6963,37.304},{11.3045,11.3041}},
        rotation=-90,
        origin={58.696,101.304})));
protected
  TransiEnt.Components.Sensors.RealGas.MassFlowSensor massflowSensor_ely(medium=medium, xiNumber=0)         annotation (Placement(transformation(
        extent={{7,6},{-7,-6}},
        rotation=180,
        origin={-2,6})));

equation

  connect(gasPortOut, massflowSensor_ely.gasPortOut) annotation (Line(
      points={{5,-94},{6,-94},{6,-4},{5,-4},{5,0}},
      color={255,255,0},
      thickness=1.5));
  connect(electrolyzer.gasPortOut, massflowSensor_ely.gasPortIn) annotation (Line(
      points={{-26,0},{-9,0}},
      color={255,255,0},
      thickness=1.5));
  connect(nulPower.y, electrolyzer.storageCompressorPowerIn) annotation (Line(points={{-81.3,-69},{-81.3,-70},{-78,-70},{-78,-13.68},{-70,-13.68}}, color={0,0,127}));
  connect(minMaxController.P_el_ely, electrolyzer.P_el_set) annotation (Line(
      points={{-54,35.2},{-53.6,36},{-53.6,21.6}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(P_el_set, minMaxController.P_el_set) annotation (Line(points={{0,108},{0,78},{-54,78},{-54,57}}, color={0,127,127}));
  connect(electrolyzer.electrolyzerPowerOut, electrolyzerPowerOut) annotation (Line(
      points={{-57.2,-12.6},{-57.2,-12},{-48,-12},{-48,-76},{104,-76}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(nulPower.y, electrolyzer.coolingPumpPowerIn) annotation (Line(points={{-81.3,-69},{-81.3,-70},{-78,-70},{-78,-52},{-16,-52},{-16,-4},{-14,-4},{-14,3.6},{-22,3.6}}, color={0,0,127}));
  connect(coolingModel_Valve.T_op, electrolyzer.temperatureOut) annotation (Line(points={{33.96,-16.42},{-14,-16.42},{-14,-20},{-80,-20},{-80,-5.76},{-57.2,-5.76}}, color={0,0,127}));
  connect(electrolyzer.heat, coolingModel_Valve.heatPortCooling) annotation (Line(points={{-26,-11.88},{24,-11.88},{24,-36.96},{36,-36.96}}, color={191,0,0}));
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
<p><span style=\"font-family: Arial;\">Electrolyzer system with L2 electrolyzer model and its accounted BoP, consisting of a cooling model with a pump. All models are found in H2Microgrid_TransiEnt and are based on existing TransiEnt models.</span></p>
<p><span style=\"font-family: Arial;\">This mode is meatnt to be used with a low-pressure storage and the compressor power is considered as nul in the electrolyzer model.</span></p>
<p><span style=\"font-family: Arial;\">The electrolyzer system also contains a power controller found in TransiEnt library.</span></p>
<p><br>Tested in the check models &quot;H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Tests.Test_PEMElectrolyzerL2_nonCompressedStorage&quot;</p>
<p><br>Simplified feed-in station model, without connection to a natural gas grid</p>
</html>"));
end SystemElectrolyzerL2_nonCompressedStorage;
