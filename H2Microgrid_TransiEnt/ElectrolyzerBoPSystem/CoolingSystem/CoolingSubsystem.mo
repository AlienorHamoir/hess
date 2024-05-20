within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.CoolingSystem;
model CoolingSubsystem
 replaceable package Coolant_Medium =
      Modelica.Media.Water.ConstantPropertyLiquidWater                                  constrainedby Modelica.Media.Interfaces.PartialMedium;
 inner Modelica.Fluid.System system(energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
      massDynamics=Modelica.Fluid.Types.Dynamics.SteadyStateInitial)                                                                                               annotation (
    Placement(visible = true, transformation(origin={-84,-78},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  VirtualFCS.Thermal.HeatSink heatSink(redeclare package Medium =
        Coolant_Medium)                                                           annotation (
    Placement(visible = true, transformation(origin={-12,-38},    extent = {{-10, -10}, {10, 10}}, rotation = 180)));
  Modelica.Blocks.Sources.RealExpression setPumpSpeed(y=subSystemCoolingControl.controlInterface)   annotation (
    Placement(visible = true, transformation(origin={16,-2},   extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  VirtualFCS.SubSystems.Cooling.SubSystemCoolingControl subSystemCoolingControl annotation (
    Placement(visible = true, transformation(origin={-32,66},    extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Fluid.Vessels.OpenTank tankCoolant(
    redeclare package Medium = Coolant_Medium,
    crossArea=0.0314,
    height=0.16,
    level_start=0.12,
    nPorts=1,
    massDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    use_HeatTransfer=true,
    portsData={Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter=0.1)},
    T_start=system.T_start)                                                                                                                                                                                                         annotation (
    Placement(visible = true, transformation(origin={-69,19},    extent = {{-11, -11}, {11, 11}}, rotation = 0)));

  VirtualFCS.Fluid.PumpElectricDC pumpElectricDC annotation (
    Placement(visible = true, transformation(origin={40,18},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Fittings.TeeJunctionVolume teeJunctionCoolantTank(redeclare package Medium =
                       Coolant_Medium, V=0.00001)                                                                          annotation (
    Placement(visible = true, transformation(origin={-36,-10},    extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium =
        Coolant_Medium)                                                                   annotation (
    Placement(visible = true, transformation(origin={108,18},    extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 120}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium =
        Coolant_Medium)                                                                   annotation (
    Placement(visible = true, transformation(origin={108,-30},    extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, 120}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation (
    Placement(visible = true, transformation(origin={58,100},   extent = {{10, -10}, {-10, 10}}, rotation = 0), iconTransformation(origin = {50, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation (
    Placement(visible = true, transformation(origin={28,100},   extent = {{10, -10}, {-10, 10}}, rotation = 0), iconTransformation(origin = {-50, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput sensors[2] annotation (
    Placement(visible = true, transformation(origin={106,-6},   extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
  Modelica.Blocks.Interfaces.RealInput controlInterface annotation (
    Placement(visible = true, transformation(origin={-106,66},    extent = {{-14, -14}, {14, 14}}, rotation = 0), iconTransformation(origin = {-110, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Machines.Pump pump annotation (Placement(transformation(extent={{-42,-104},{-22,-84}})));
  Buildings.Experimental.DHC.EnergyTransferStations.BaseClasses.Pump_m_flow pump_m_flow annotation (Placement(transformation(extent={{38,-108},{58,-88}})));
  Buildings.Experimental.DHC.Loads.HotWater.BaseClasses.HeatExchangerPumpController conPum annotation (Placement(transformation(extent={{0,-106},{20,-86}})));
equation
  connect(heatSink.port_b,teeJunctionCoolantTank. port_1) annotation (
    Line(points={{-22,-31},{-36,-31},{-36,-20}},        color = {0, 0, 255}, thickness = 1.5));
  connect(tankCoolant.ports[1],teeJunctionCoolantTank. port_3) annotation (
    Line(points={{-69,8},{-69,-10},{-46,-10}},        color = {0, 0, 255}, thickness = 1.5));
  connect(teeJunctionCoolantTank.port_2,pumpElectricDC. Input) annotation (
    Line(points={{-36,0},{-36,18},{31,18}},         color = {0, 0, 255}, thickness = 1.5));
  connect(pin_n,pumpElectricDC. pin_n) annotation (
    Line(points={{28,100},{28,32},{37,32},{37,26}},         color = {0, 0, 255}));
  connect(setPumpSpeed.y,pumpElectricDC. contol_input) annotation (
    Line(points={{27,-2},{37,-2},{37,10}},       color = {0, 0, 127}));
  connect(pumpElectricDC.sensors,sensors)  annotation (
    Line(points={{43,10},{43,-6},{106,-6}},      color = {0, 0, 127}, thickness = 0.5));
  connect(port_a,heatSink. port_a) annotation (
    Line(points={{108,-30},{54,-30},{54,-31},{-2,-31}},
                                            color = {255, 0, 0}, thickness = 1.5));
  connect(pumpElectricDC.Output,port_b)  annotation (
    Line(points={{49,18},{108,18}},      color = {0, 0, 255}, thickness = 1.5));
  connect(controlInterface,subSystemCoolingControl. sensorInterface) annotation (
    Line(points={{-106,66},{-54,66}},      color = {0, 0, 127}));
  connect(pin_p, pin_p)
    annotation (Line(points={{58,100},{58,100}}, color={0,0,255}));
  connect(pumpElectricDC.pin_p, pin_p) annotation (Line(points={{43,26},{43,86},
          {58,86},{58,100}}, color={0,0,255}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={102,44,145},
          fillColor={18,69,127},
          fillPattern=FillPattern.Solid),                                                                                                                                                               Text(origin={-48,70},    textColor={255,255,
              255},                                                                                                                                                                                                        extent = {{-22, 12}, {112, -142}},
          textString="Cool")}), Diagram(coordinateSystem(preserveAspectRatio=false),
        graphics={      Text(origin={45,30},    extent = {{-19, 4}, {15, -2}}, textString = "Pump")}));
end CoolingSubsystem;
