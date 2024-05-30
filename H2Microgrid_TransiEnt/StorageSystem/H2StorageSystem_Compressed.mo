within H2Microgrid_TransiEnt.StorageSystem;
model H2StorageSystem_Compressed
protected
  TransiEnt.Components.Gas.Compressor.Controller.ControllerCompressor_dp controlCompressorAfterELY(
    p_paramBefore=true,
    p_paramAfter=false,
    p_beforeCompParam=20e5,
    p_afterCompParam=30e5)                                                                                              annotation (Placement(transformation(extent={{-24,20},{-4,40}})));
public
  TransiEnt.Storage.Gas.GasStorage_constXi_L2 storage(
    start_pressure=start_pressure,
    includeHeatTransfer=includeHeatTransfer,
    medium=medium,
    V_geo=V_geo,
    redeclare model HeatTransfer = TransiEnt.Storage.Gas.Base.ConstantHTOuterTemperature_L2 (alpha_nom=alpha_nom),
    m_gas_start=m_start,
    p_gas_start=p_start,
    T_gas_start=T_start,
    p_max=p_maxHigh,
    eta_ely=eta_n,
    Cspec_demAndRev_el=Cspec_demAndRev_el_other) annotation (Placement(transformation(
        extent={{14,-15},{-14,15}},
        rotation=90,
        origin={22,-25})));

  replaceable parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3 "Hydrogen model to be used" annotation (Dialog(tab="General", group="General"));

  parameter Boolean start_pressure=true "true if a start pressure is defined, false if a start mass is defined" annotation (Dialog(tab="Storage"));
  parameter Boolean includeHeatTransfer=true "false for neglecting heat transfer" annotation (Dialog(tab="Storage"));
  parameter Modelica.Units.SI.Volume V_geo=0.1 "Geometric volume of storage" annotation (Dialog(tab="Storage"));
  parameter Modelica.Units.SI.Height height=3.779*V_geo^(1/3) "Height of storage" annotation (Dialog(tab="Storage"));
  parameter Modelica.Units.SI.CoefficientOfHeatTransfer alpha_nom=4 "Heat transfer coefficient inside the storage cylinder" annotation (Dialog(tab="Storage"));
  parameter Modelica.Units.SI.Mass m_start=1 "Stored gas mass at t=0" annotation (Dialog(tab="Storage"));
  parameter Modelica.Units.SI.Pressure p_start=simCenter.p_amb_const + simCenter.p_eff_2  "Pressure in storage at t=0" annotation (Dialog(tab="Storage")); // p_start = 17 bar
  parameter Modelica.Units.SI.ThermodynamicTemperature T_start=283.15 "Temperature of gas in storage at t=0" annotation (Dialog(tab="Storage"));
  parameter Modelica.Units.SI.Pressure p_out = 20e5;
  parameter Modelica.Units.SI.Pressure p_maxLow=p_maxHigh - 1e5 "Lower limit of the maximum pressure in storage" annotation (Dialog(tab="Storage", group="Control"));
  parameter Modelica.Units.SI.Pressure p_maxHigh=30e5 "Upper limit of the maximum pressure in storage" annotation (Dialog(tab="Storage", group="Control"));
  parameter Modelica.Units.SI.Pressure p_minLow_constantDemand=0.5e3 "storage can be emptied via 'm_flow_hydrogenDemand_constant' up to 'p_minLow_constantDemand'" annotation (Dialog(tab="Storage", group="Control"));
  parameter Modelica.Units.SI.MassFlowRate m_flow_hydrogenDemand_constant=0 "constant hydrogen demand if hydrogen is available" annotation (Dialog(tab="Storage", group="Control"));

 parameter Modelica.Units.SI.Temperature T_out=283.15 "Hydrogen output temperature" annotation(Dialog(group="Fundamental Definitions"));

  parameter Modelica.Units.SI.Pressure dp_Low=1e3 "if valve open and (p_in-p_out) <= dp_Low, close valve" annotation (Dialog(tab="Storage", group="Control"));
  parameter Modelica.Units.SI.Pressure dp_High=1e5 "if valve closed and (p_in-p_out) >= dp_High, open valve" annotation (Dialog(tab="Storage", group="Control"));

  parameter Boolean StoreAllHydrogen=false "|Storage|Control|All Hydrogen is stored before beeing fed in" annotation (Dialog(tab="Storage", group="Control"));

  parameter Modelica.Units.SI.Efficiency eta_mech_compressor(
    min=0,
    max=1)=0.8 "Sets a with increasing input power linear degrading efficiency coefficient (min = 0, max = 1)" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.Efficiency eta_el_compressor(
    min=0,
    max=1)=0.8 "Sets a with increasing input power linear degrading efficiency coefficient (min = 0, max = 1)" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.Power P_el_n_compressor=200;

  parameter Modelica.Units.SI.Efficiency eta_n(
    min=0,
    max=1)=0.75 "Nominal efficency refering to the GCV (min = 0, max = 1)" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter TransiEnt.Basics.Units.MonetaryUnitPerEnergy Cspec_demAndRev_el_other=1;

  TransiEnt.Basics.Interfaces.Gas.RealGasPortIn H2PortIn(Medium=medium) "inlet flow" annotation (Placement(transformation(extent={{-12,-114},{14,-88}}), iconTransformation(extent={{-12,-114},{14,-88}})));
  TransiEnt.Basics.Interfaces.Gas.RealGasPortOut H2PortOut(Medium=medium) "outlet flow" annotation (Placement(transformation(extent={{-14,84},{14,112}}), iconTransformation(extent={{-14,84},{14,112}})));
  TransiEnt.Basics.Interfaces.General.PressureOut pressureTank annotation (Placement(transformation(extent={{96,50},{116,70}})));
  StorageSystem.TankSOC tankSOC annotation (Placement(transformation(extent={{42,76},{62,96}})));
  Modelica.Blocks.Interfaces.RealOutput socTank annotation (Placement(transformation(extent={{96,76},{116,96}})));

// public
//   inner Summary summary(
//     storage(
//       mediumModel=storage.summary.gasBulk.mediumModel,
//       xi=storage.summary.gasBulk.xi,
//       x=storage.summary.gasBulk.x,
//       mass=storage.summary.gasBulk.mass,
//       T=storage.summary.gasBulk.T,
//       p=storage.summary.gasBulk.p,
//       h=storage.summary.gasBulk.h,
//       rho=storage.summary.gasBulk.rho),
//     gasPortStorageIn(
//       mediumModel=storage.summary.gasPortIn.mediumModel,
//       xi=storage.summary.gasPortIn.xi,
//       x=storage.summary.gasPortIn.x,
//       m_flow=storage.summary.gasPortIn.m_flow,
//       T=storage.summary.gasPortIn.T,
//       p=storage.summary.gasPortIn.p,
//       h=storage.summary.gasPortIn.h,
//       rho=storage.summary.gasPortIn.rho),
//     gasPortStorageOut(
//       mediumModel=storage.summary.gasPortOut.mediumModel,
//       xi=storage.summary.gasPortOut.xi,
//       x=storage.summary.gasPortOut.x,
//       m_flow=storage.summary.gasPortOut.m_flow,
//       T=storage.summary.gasPortOut.T,
//       p=storage.summary.gasPortOut.p,
//       h=storage.summary.gasPortOut.h,
//       rho=storage.summary.gasPortOut.rho),
//     costsStorage(
//       costs=storage.summary.costs.costs,
//       investCosts=storage.summary.costs.investCosts,
//       investCostsStartGas=storage.summary.costs.investCostsStartGas,
//       demandCosts=storage.summary.costs.demandCosts,
//       oMCosts=storage.summary.costs.oMCosts,
//       otherCosts=storage.summary.costs.otherCosts,
//       revenues=storage.summary.costs.revenues));

  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{64,-94},{84,-74}})));
  TransiEnt.Components.Electrical.Machines.MotorComplex
                                   motorComplex(eta=0.97) annotation (Placement(transformation(extent={{-6,-74},{14,-54}})));
  TransiEnt.Components.Boundaries.Electrical.ComplexPower.SlackBoundary
                                                   slackBoundary annotation (Placement(transformation(extent={{24,-74},{44,-54}})));
  TransiEnt.Components.Gas.Compressor.CompressorRealGasIsothermal_L1_simple compressorRealGasIsothermal_L1_simple(
    medium=medium,
    presetVariableType="dp",
    useMechPowerPort=false,
    use_Delta_p_input=true) annotation (Placement(transformation(extent={{-74,-36},{-54,-16}})));
  TransiEnt.Components.Sensors.RealGas.MassFlowSensor massFlowSensorToFC(medium=medium) annotation (Placement(transformation(extent={{58,-26},{78,-6}})));
  TransiEnt.Components.Sensors.RealGas.MassFlowSensor massFlowSensorFromElectrolyzer(medium=medium) annotation (Placement(transformation(extent={{-32,-26},{-12,-6}})));
equation
  connect(tankSOC.tankSOC, socTank) annotation (Line(points={{62.6,86},{106,86}}, color={0,0,127}));
  connect(motorComplex.mpp, compressorRealGasIsothermal_L1_simple.mpp) annotation (Line(points={{-6,-64},{-64,-64},{-64,-36}},  color={95,95,95}));
  connect(slackBoundary.epp, motorComplex.epp) annotation (Line(
      points={{24,-64},{24,-64.1},{14.1,-64.1}},
      color={28,108,200},
      thickness=0.5));
  connect(controlCompressorAfterELY.Delta_p, compressorRealGasIsothermal_L1_simple.dp_in) annotation (Line(points={{-14,19},{-14,-8},{-56,-8},{-56,-15}},
                                                                                                                                       color={0,0,127}));
  connect(storage.gasPortOut, massFlowSensorToFC.gasPortIn) annotation (Line(
      points={{31.45,-25},{31.45,-26},{58,-26}},
      color={255,255,0},
      thickness=1.5));
  connect(massFlowSensorToFC.gasPortOut, H2PortOut) annotation (Line(
      points={{78,-26},{78,-28},{88,-28},{88,98},{0,98}},
      color={255,255,0},
      thickness=1.5));
  connect(H2PortIn, compressorRealGasIsothermal_L1_simple.gasPortIn) annotation (Line(
      points={{1,-101},{-44,-101},{-44,-102},{-88,-102},{-88,-26},{-74,-26}},
      color={255,255,0},
      thickness=1.5));
  connect(compressorRealGasIsothermal_L1_simple.gasPortOut, massFlowSensorFromElectrolyzer.gasPortIn) annotation (Line(
      points={{-54,-26},{-32,-26}},
      color={255,255,0},
      thickness=1.5));
  connect(massFlowSensorFromElectrolyzer.gasPortOut, storage.gasPortIn) annotation (Line(
      points={{-12,-26},{2,-26},{2,-25},{14.65,-25}},
      color={255,255,0},
      thickness=1.5));
  connect(storage.p_gas, controlCompressorAfterELY.p_afterCompIn) annotation (Line(points={{22,-18},{22,30},{-4,30}}, color={0,0,127}));
  connect(storage.p_gas, tankSOC.currentTankPressure) annotation (Line(points={{22,-18},{22,85.8},{41.8,85.8}}, color={0,0,127}));
  connect(storage.p_gas, pressureTank) annotation (Line(points={{22,-18},{22,60},{106,60}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-60,-62},{60,-82}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid),
        Line(points={{-60,68},{-60,-72}}, color={28,108,200}),
        Ellipse(
          extent={{-60,78},{60,58}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid),
        Line(points={{60,68},{60,-72}}, color={28,108,200})}), Diagram(coordinateSystem(preserveAspectRatio=false)));
end H2StorageSystem_Compressed;
