within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.Systems;
model SystemElectrolyzerL1_nonCompressedStorage

  extends TransiEnt.Producer.Gas.Electrolyzer.Base.PartialFeedInStation(
    gasPortOut(Medium=medium) annotation(Placement(
      transformation(
        origin={400, 400},
        extent={{-10, -10}, {10, 10}},
        rotation=0))),
    break m_flow_feedIn);

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_coolant=simCenter.fluid1;

  parameter Modelica.Units.SI.ActivePower P_el_n=5e3 "nominal power of electrolyser" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.ActivePower P_el_max=1.0*P_el_n "Maximum power of electrolyzer" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.ActivePower P_el_min=0.1*P_el_n "Minimal power of electrolyzer" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.MassFlowRate m_flow_start=0.0 "Sets initial value for m_flow from a buffer" annotation (Dialog(tab="General", group="Initialization"));
  // parameter SI.Temperature T_Init=283.15 "Sets initial value for T" annotation (Dialog(tab="General", group="Initialization"));
//  parameter Real eta_inv_n=0.956 "nominal efficiency of the inverter" annotation (Dialog(group="Fundamental Definitions")); // we neglect all losses relative to grid constraints, including inverter efficiency
  parameter Modelica.Units.SI.Power P_el_pump = 285.7 "W, pump el power consumption";

  parameter Modelica.Units.SI.Efficiency eta_n(
    min=0,
    max=1)=0.755 "Nominal efficiency refering to the GCV (min = 0, max = 1)" annotation (Dialog(group="Fundamental Definitions"));
  parameter Modelica.Units.SI.Efficiency eta_scale(
    min=0,
    max=1)=0.01 "Sets a with increasing input power linear degrading efficiency coefficient (min=0,max=1)" annotation (Dialog(group="Fundamental Definitions"));

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
  PEMElectrolyzer_L1 electrolyzer1(
    useHeatPort=false,
    integrateH2Flow=true,
    integrateElPower=true,
    calculateCost=false,
    T_out=313.15,
    T_amb(displayUnit="K")=273.15+23,
    specificWaterConsumption=11,
    useLeakageMassFlow=false,
    redeclare model Dynamics = TransiEnt.Producer.Gas.Electrolyzer.Base.ElectrolyzerDynamics1stOrder,
    redeclare model Charline = TransiEnt.Producer.Gas.Electrolyzer.Base.ElectrolyzerEfficiencyCharlineSilyzer200)
                              annotation (Placement(transformation(extent={{-62,-16},{-30,16}})));
  Modelica.Blocks.Sources.Constant nulCompressorPower(k=0) annotation (Placement(transformation(extent={{-92,-58},{-78,-44}})));
  TransiEnt.Producer.Gas.Electrolyzer.Controller.MinMaxController minMaxController(
    P_el_n=P_el_n,
    P_el_max=P_el_max + P_el_pump,
    P_el_min=P_el_min + P_el_pump) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={-52,52})));
  Modelica.Blocks.Interfaces.RealOutput electrolyzerPowerOut "Output signal connector" annotation (Placement(transformation(extent={{92,-70},{112,-50}})));
protected
  TransiEnt.Components.Sensors.RealGas.MassFlowSensor massflowSensor_ely(medium=medium, xiNumber=0)         annotation (Placement(transformation(
        extent={{7,6},{-7,-6}},
        rotation=180,
        origin={-2,6})));

equation
  connect(massflowSensor_ely.m_flow, H2massFlowRateOutElectrolyzer) annotation (Line(points={{5.7,6},{88,6},{88,65},{105,65}},
                                                                                                                        color={0,0,127}));

  connect(gasPortOut, massflowSensor_ely.gasPortOut) annotation (Line(
      points={{5,-94},{6,-94},{6,0},{5,0}},
      color={255,255,0},
      thickness=1.5));
  connect(electrolyzer1.epp, epp) annotation (Line(
      points={{-62,0},{-100,0}},
      color={0,135,135},
      thickness=0.5));
  connect(electrolyzer1.gasPortOut, massflowSensor_ely.gasPortIn) annotation (Line(
      points={{-30,0},{-9,0}},
      color={255,255,0},
      thickness=1.5));
  connect(nulCompressorPower.y, electrolyzer1.storageCompressorPowerIn) annotation (Line(points={{-77.3,-51},{-77.3,-52},{-72,-52},{-72,-12},{-65.2,-12},{-65.2,-12.16}}, color={0,0,127}));
  connect(nulCompressorPower.y, electrolyzer1.coolingPumpPowerIn) annotation (Line(points={{-77.3,-51},{-77.3,-52},{-62,-52},{-62,-34},{-16,-34},{-16,3.2},{-26.8,3.2}}, color={0,0,127}));
  connect(minMaxController.P_el_ely, electrolyzer1.P_el_set) annotation (Line(
      points={{-52,41.2},{-52.4,40},{-52.4,19.2}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(P_el_set, minMaxController.P_el_set) annotation (Line(points={{0,108},{0,74},{-52,74},{-52,63}}, color={0,127,127}));
  connect(electrolyzer1.electrolyzerPowerOut, electrolyzerPowerOut) annotation (Line(
      points={{-53.04,-7.68},{8,-7.68},{8,-8},{66,-8},{66,-60},{102,-60}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid), Text(
          extent={{-76,38},{76,-42}},
          textColor={255,255,255},
          textStyle={TextStyle.Bold},
          textString="Electrolyzer
L1")}),                                     Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p><span style=\"font-family: Arial;\">Electrolyzer system with L1 simplified electrolyzer model. No BoP is accounted for. All models are found in H2Microgrid_TransiEnt and are based on existing TransiEnt models.</span></p>
<p><span style=\"font-family: Arial;\">The electrolyzer system also contains a power controller found in TransiEnt library.</span></p>
<p><br>Tested in the check models &quot;H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Tests.Test_PEMElectrolyzerL1_nonCompressedStorage&quot;</p>
<p><br>Simplified feed-in station model, without connection to a natural gas grid</p>
</html>"));
end SystemElectrolyzerL1_nonCompressedStorage;
