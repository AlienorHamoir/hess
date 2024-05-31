within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.CoolingSystem.HeatPortCooling;
model Cooling_VirtualFCS

    replaceable package medium_coolant = Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby Modelica.Media.Interfaces.PartialMedium;

  VirtualFCS.Thermal.HeatSink heatSink(redeclare package Medium = medium_coolant) annotation (
    Placement(visible = true, transformation(origin={31,-61},     extent={{-7,-7},{7,7}},          rotation=270)));
  Modelica.Blocks.Sources.RealExpression setPumpSpeed(y=subSystemCoolingControl.controlInterface)   annotation (
    Placement(visible = true, transformation(origin={-4,20},   extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  VirtualFCS.SubSystems.Cooling.SubSystemCoolingControl subSystemCoolingControl annotation (
    Placement(visible = true, transformation(origin={0,64},      extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Fluid.Vessels.OpenTank tankCoolant(
    redeclare package Medium = medium_coolant,
    crossArea=0.0314,
    height=0.16,
    level_start=0.12,
    nPorts=1,
    massDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    use_HeatTransfer=true,
    portsData={Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter=0.1)},
    T_start=system.T_start)                                                                                                                                                                                                         annotation (
    Placement(visible = true, transformation(origin={-43,3},     extent = {{-11, -11}, {11, 11}}, rotation = 0)));
  VirtualFCS.Fluid.PumpElectricDC pumpElectricDC annotation (
    Placement(visible = true, transformation(origin={26,-10},   extent = {{-10, -10}, {10, 10}}, rotation=-90)));
  Modelica.Fluid.Fittings.TeeJunctionVolume teeJunctionCoolantTank(redeclare package Medium = medium_coolant, V=0.00001)   annotation (
    Placement(visible = true, transformation(origin={-22,-20},    extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Fluid.Pipes.DynamicPipe pipeCoolant(
    redeclare package Medium = medium_coolant,
    diameter=0.003,
    length=1,
    modelStructure=Modelica.Fluid.Types.ModelStructure.a_vb,
    nNodes=1,
    nParallel=500,
    p_a_start=102502,
    use_HeatTransfer=true)                                                                                                                                                                                                         annotation (
    Placement(visible = true, transformation(origin={27,-39},   extent={{-5,5},{5,-5}},          rotation=-90)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=1777)   annotation (
    Placement(visible = true, transformation(origin={-50,-40},  extent={{-8,-8},{8,8}},          rotation=0)));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor1(
    C=0.1,
    T(fixed=false, start=293.15),
    der_T(fixed=false))                                                                                                                                 annotation (
    Placement(visible = true, transformation(                    extent={{-82,-40},{-66,-24}},    rotation = 0)));
  Modelica.Electrical.Analog.Sensors.PowerSensor power_coolingPump annotation (Placement(visible=true, transformation(extent={{10,10},{-10,-10}},rotation=180,
        origin={80,-52})));
  VirtualFCS.Electrochemical.Battery.BatterySystem
                                        batterySystem(
    C_bat_pack=400,
    SOC_init=0.9,
    V_max_bat_pack=27,
    V_min_bat_pack=23,
    V_nom_bat_pack=25,
    m_bat_pack=1)                                                                                                                                                     annotation (
    Placement(visible = true, transformation(origin={82,-16},   extent = {{-10, -10}, {10, 10}}, rotation=90)));
  Modelica.Blocks.Interfaces.RealOutput sensors[2] annotation (
    Placement(visible = true, transformation(origin={106,46},   extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin={98,-12},      extent = {{-10, -10}, {10, 10}}, rotation=0)));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPortElectrolyzer annotation (Placement(transformation(extent={{-110,-50},{-90,-30}})));
  Modelica.Blocks.Interfaces.RealInput temperatureElectrolyzer annotation (Placement(transformation(extent={{-128,44},{-88,84}})));
  Modelica.Blocks.Interfaces.RealOutput powerCoolingCompressor "Instantaneous power as output signal" annotation (Placement(transformation(extent={{92,-102},{124,-70}}), iconTransformation(extent={{92,-102},{124,-70}})));
equation
  connect(heatSink.port_b,teeJunctionCoolantTank. port_1) annotation (
    Line(points={{26.1,-68},{26,-68},{26,-74},{-22,-74},{-22,-30}},
                                                        color = {0, 0, 255}, thickness = 1.5));
  connect(tankCoolant.ports[1],teeJunctionCoolantTank. port_3) annotation (
    Line(points={{-43,-8},{-44,-8},{-44,-20},{-32,-20}},
                                                      color = {0, 0, 255}, thickness = 1.5));
  connect(teeJunctionCoolantTank.port_2,pumpElectricDC. Input) annotation (
    Line(points={{-22,-10},{-22,8},{26,8},{26,-1}}, color = {0, 0, 255}, thickness = 1.5));
  connect(setPumpSpeed.y,pumpElectricDC. contol_input) annotation (
    Line(points={{7,20},{12,20},{12,-7},{18,-7}},color = {0, 0, 127}));
  connect(pumpElectricDC.sensors,sensors)  annotation (
    Line(points={{18,-13},{18,-12},{8,-12},{8,12},{52,12},{52,46},{106,46}},
                                                 color = {0, 0, 127}, thickness = 0.5));
  connect(thermalConductor.port_b,pipeCoolant. heatPorts[1]) annotation (
    Line(points={{-42,-40},{6,-40},{6,-39.05},{24.8,-39.05}},
                                                            color = {191, 0, 0}));
  connect(thermalConductor.port_a,heatCapacitor1. port) annotation (Line(points={{-58,-40},{-74,-40}},
                                                                                                   color={191,0,0}));
  connect(pumpElectricDC.Output,pipeCoolant. port_a) annotation (Line(points={{26,-19},{27,-19},{27,-34}},
                                                                                                        color={0,127,255}));
  connect(heatSink.port_a,pipeCoolant. port_b) annotation (Line(points={{26.1,-54},{27,-54},{27,-44}},
                                                                                                    color={0,127,255}));
  connect(heatCapacitor1.port, heatPortElectrolyzer) annotation (Line(points={{-74,-40},{-100,-40}}, color={191,0,0}));
  connect(subSystemCoolingControl.sensorInterface, temperatureElectrolyzer) annotation (Line(points={{-22,64},{-108,64}}, color={0,0,127}));
  connect(power_coolingPump.power, powerCoolingCompressor) annotation (Line(points={{70,-63},{70,-86},{108,-86}}, color={0,0,127}));
  connect(power_coolingPump.pv, power_coolingPump.pc) annotation (Line(points={{80,-42},{66,-42},{66,-52},{70,-52}}, color={0,0,255}));
  connect(power_coolingPump.nv, power_coolingPump.nc) annotation (Line(points={{80,-62},{80,-66},{94,-66},{94,-52},{90,-52}}, color={0,0,255}));
  connect(pumpElectricDC.pin_p, batterySystem.pin_n) annotation (Line(points={{34,-13},{34,-14},{66,-14},{66,-20.4},{72.4,-20.4}}, color={0,0,255}));
  connect(pumpElectricDC.pin_n, batterySystem.pin_p) annotation (Line(points={{34,-7},{34,-8},{66,-8},{66,-11.6},{72.4,-11.6}}, color={0,0,255}));
  connect(power_coolingPump.pc, pumpElectricDC.pin_p) annotation (Line(points={{70,-52},{66,-52},{66,-13},{34,-13}}, color={0,0,255}));
  connect(power_coolingPump.nv, pumpElectricDC.pin_n) annotation (Line(points={{80,-62},{80,-66},{94,-66},{94,-30},{96,-30},{96,-2},{46,-2},{46,-7},{34,-7}}, color={0,0,255}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end Cooling_VirtualFCS;
