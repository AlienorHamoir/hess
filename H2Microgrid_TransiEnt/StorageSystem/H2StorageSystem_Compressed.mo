within H2Microgrid_TransiEnt.StorageSystem;
model H2StorageSystem_Compressed "High-pressure compressed storage"
protected
  TransiEnt.Components.Gas.Compressor.Controller.ControllerCompressor_dp controlCompressor(
    p_paramBefore=true,
    p_paramAfter=false,
    p_beforeCompParam=p_out,
    p_afterCompParam=p_maxHigh) annotation (Placement(transformation(extent={{-74,46},{-54,66}})));
public
  TransiEnt.Storage.Gas.GasStorage_constXi_L2 H2storage(
    calculateCost=false,
    start_pressure=start_pressure,
    includeHeatTransfer=false,
    medium=medium,
    V_geo=V_geo,
    redeclare model HeatTransfer = TransiEnt.Storage.Gas.Base.ConstantHTOuterTemperature_L2 (alpha_nom=alpha_nom),
    m_gas_start=m_start,
    p_gas_start=p_start,
    T_gas_start=T_start,
    p_max=p_maxHigh,
    eta_ely=eta_n) annotation (Placement(transformation(
        extent={{15.5,-17.5},{-15.5,17.5}},
        rotation=90,
        origin={-3.5,0.5})));

  replaceable parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3 "Hydrogen model to be used" annotation (Dialog(tab="General", group="General"));

  parameter Boolean start_pressure=true "true if a start pressure is defined, false if a start mass is defined" annotation (Dialog(tab="General", group="Storage"));
  parameter Boolean includeHeatTransfer=false "false for neglecting heat transfer" annotation (Dialog(tab="General", group="Storage"));
  parameter Modelica.Units.SI.Volume V_geo=1.6 "m3, Geometric volume of storage" annotation (Dialog(tab="General", group="Storage"));
  parameter Modelica.Units.SI.Height height=3.779*V_geo^(1/3) "Height of storage" annotation (Dialog(tab="General", group="Storage"));
  parameter Modelica.Units.SI.CoefficientOfHeatTransfer alpha_nom=4 "Heat transfer coefficient inside the storage cylinder" annotation (Dialog(tab="General", group="Storage"));
  parameter Modelica.Units.SI.Mass m_start=1 "Stored gas mass at t=0" annotation (Dialog(tab="General", group="Storage"));
  parameter Modelica.Units.SI.Pressure p_start=15e5   "Pressure in storage at t=0" annotation (Dialog(tab="General", group="Storage"));
  parameter Modelica.Units.SI.ThermodynamicTemperature T_start=T_out "Temperature of gas in storage at t=0" annotation (Dialog(tab="General", group="Storage"));
  parameter Modelica.Units.SI.Pressure p_out = 30e5 "Output pressure of hydrogen from electrolyzer" annotation (Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.Pressure p_maxLow=p_maxHigh - 1e5 "Lower limit of the maximum pressure in storage" annotation (Dialog(tab="General", group="Control"));
  parameter Modelica.Units.SI.Pressure p_maxHigh=350e5 "Upper limit of the maximum pressure in storage" annotation (Dialog(tab="General", group="Control"));


  parameter Modelica.Units.SI.Temperature T_out=273.15+40 "Hydrogen output temperature" annotation(Dialog(tab="General", group="Electrolyzer"));
  parameter Modelica.Units.SI.Efficiency eta_mech_compressor(
    min=0,
    max=1)=0.9 "Compressor mechanical efficiency coefficient (min = 0, max = 1)" annotation (Dialog(tab="General", group="Compressor"));
  parameter Modelica.Units.SI.Efficiency eta_el_compressor(
    min=0,
    max=1)=0.9 "Compressor motor electrical efficiency coefficient (min = 0, max = 1)" annotation (Dialog(tab="General", group="Compressor"));
  parameter Modelica.Units.SI.Power P_el_n_compressor=62.72 "W, compressor nominal electrical power" annotation (Dialog(tab="General", group="Compressor"));
  parameter Modelica.Units.SI.Efficiency eta_n(
    min=0,
    max=1)=0.743 "Nominal efficency refering to the GCV (min = 0, max = 1)" annotation (Dialog(tab="General", group="Electrolyzer"));

  TransiEnt.Basics.Interfaces.Gas.RealGasPortIn H2PortIn(Medium=medium) "inlet flow" annotation (Placement(transformation(extent={{-12,-114},{14,-88}}), iconTransformation(extent={{-12,-114},{14,-88}})));
  TransiEnt.Basics.Interfaces.Gas.RealGasPortOut H2PortOut(Medium=medium) "outlet flow" annotation (Placement(transformation(extent={{-14,84},{14,112}}), iconTransformation(extent={{-14,84},{14,112}})));
  TransiEnt.Basics.Interfaces.General.PressureOut pressureTank annotation (Placement(transformation(extent={{96,30},{116,50}})));
  StorageSystem.TankSOC tankSOC(maxPressure=p_maxHigh)
                                annotation (Placement(transformation(extent={{42,68},{62,88}})));
  Modelica.Blocks.Interfaces.RealOutput socTank annotation (Placement(transformation(extent={{96,68},{116,88}})));

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

  TransiEnt.Components.Electrical.Machines.MotorComplex
                                   motorComplex(cosphi=1, eta=eta_el_compressor)
                                                          annotation (Placement(transformation(extent={{-46,-60},{-26,-40}})));
  TransiEnt.Components.Boundaries.Electrical.ComplexPower.SlackBoundary
                                                   slackBoundary annotation (Placement(transformation(extent={{32,-60},{52,-40}})));
  TransiEnt.Components.Sensors.ElectricPowerComplex electricPowerComplex annotation (Placement(transformation(extent={{-6,-60},{14,-40}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_comp "Compressor active Power" annotation (Placement(transformation(extent={{94,-46},{114,-26}})));
  ValveAndCompressor_dp ValveAndCompressor(
    medium=medium,
    p_startSplit=simCenter.p_amb_const,
    p_startJunction=simCenter.p_amb_const,
    T_startSplit=T_start,
    T_startJunction=T_start,
    redeclare model Compressor = TransiEnt.Components.Gas.Compressor.CompressorRealGasIsentropicEff_L1_simple) annotation (Placement(transformation(extent={{-78,-14},{-48,16}})));
  inner TransiEnt.ModelStatistics modelStatistics annotation (Placement(transformation(extent={{30,-96},{50,-76}})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{66,-98},{86,-78}})));
equation
  connect(tankSOC.tankSOC, socTank) annotation (Line(points={{62.6,78},{106,78}}, color={0,0,127}));
  connect(H2storage.p_gas, controlCompressor.p_afterCompIn) annotation (Line(points={{-3.5,8.25},{-3.5,56},{-54,56}},
                                                                                                                 color={0,0,127}));
  connect(H2storage.p_gas, tankSOC.currentTankPressure) annotation (Line(points={{-3.5,8.25},{-3.5,56},{32,56},{32,78},{36,78},{36,77.8},{41.8,77.8}},
                                                                                                                  color={0,0,127}));
  connect(H2storage.p_gas, pressureTank) annotation (Line(points={{-3.5,8.25},{-3.5,56},{32,56},{32,40},{106,40}},
                                                                                              color={0,0,127}));
  connect(motorComplex.epp, electricPowerComplex.epp_IN) annotation (Line(
      points={{-25.9,-50.1},{-24,-50},{-5.2,-50}},
      color={28,108,200},
      thickness=0.5));
  connect(electricPowerComplex.epp_OUT, slackBoundary.epp) annotation (Line(
      points={{13.4,-50},{32,-50}},
      color={28,108,200},
      thickness=0.5));
  connect(H2storage.gasPortOut, H2PortOut) annotation (Line(
      points={{7.525,0.5},{7.525,0},{12,0},{12,80},{0,80},{0,98}},
      color={255,255,0},
      thickness=1.5));
  connect(pressureTank, pressureTank) annotation (Line(points={{106,40},{106,40}}, color={0,0,127}));
  connect(ValveAndCompressor.gasPortIn, H2PortIn) annotation (Line(
      points={{-78,1},{-94,1},{-94,-84},{1,-84},{1,-101}},
      color={255,255,0},
      thickness=1.5));
  connect(ValveAndCompressor.gasPortOut, H2storage.gasPortIn) annotation (Line(
      points={{-48,1},{-32,1},{-32,0.5},{-12.075,0.5}},
      color={255,255,0},
      thickness=1.5));
  connect(ValveAndCompressor.dp_desired, controlCompressor.Delta_p) annotation (Line(points={{-63,16},{-64,16},{-64,45}}, color={0,0,127}));
  connect(ValveAndCompressor.mpp2, motorComplex.mpp) annotation (Line(points={{-53.1,-13.4},{-53.1,-50},{-46,-50}}, color={95,95,95}));
  connect(electricPowerComplex.P, P_comp) annotation (Line(
      points={{-1,-41.4},{-2,-41.4},{-2,-36},{104,-36}},
      color={0,135,135},
      pattern=LinePattern.Dash));
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
        Line(points={{60,68},{60,-72}}, color={28,108,200})}), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<h4>1. Purpose of model</h4>
<p>High-pressure storage system fo HESS applications (storing hydrogen between an electrolyzer and a fuel cell). </p>
<p>Contains the model of a simple gas storage volume for constant composition, a compressor and valve model, as well as a compressor and valve controller. Compressor power is computed via the motor.</p>
<p>Tank state of charge (SOC) is also computed through a small model.</p>
<h4>2. Level of detail, physical effects considered, and physical insight</h4>
<h4>3. Limits of validity </h4>
<h4>4. Interfaces</h4>
<h4>5. Nomenclature</h4>
<p>(no remarks)</p>
<h4>6. Governing Equations</h4>
<h4>7. Remarks for Usage</h4>
<h4>8. Validation</h4>
<p>Tested in &quot;H2Microgrid_TransiEnt.StorageSystem.TestStorage_CompSystem&quot;</p>
<h4>9. References</h4>
<p>Source: initial storage model in feed-in station from TransiEnt library &quot;TransiEnt.Producer.Gas.Electrolyzer.Systems.FeedInStation_Storage&quot;</p>
<h4>10. Version History</h4>
<p>Model created for HESS applications by Ali&eacute;nor Hamoir in June 2024</p>
</html>"));
end H2StorageSystem_Compressed;
