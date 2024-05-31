within H2Microgrid_TransiEnt.StorageSystem;
model H2StorageSystem_nonCompressed
public
  TransiEnt.Storage.Gas.GasStorage_constXi_L2 storage(
    calculateCost=false,
    start_pressure=start_pressure,
    includeHeatTransfer=includeHeatTransfer,
    medium=medium,
    V_geo=V_geo,
    redeclare model HeatTransfer = TransiEnt.Storage.Gas.Base.ConstantHTOuterTemperature_L2 (alpha_nom=alpha_nom),
    m_gas_start=m_start,
    p_gas_start=p_start,
    T_gas_start=T_start,
    p_max=p_maxHigh,
    eta_ely=eta_n) annotation (Placement(transformation(
        extent={{14,-15},{-14,15}},
        rotation=90,
        origin={0,1})));

  replaceable parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3 "Hydrogen model to be used" annotation (Dialog(tab="General", group="General"));

  parameter Boolean start_pressure=true "true if a start pressure is defined, false if a start mass is defined" annotation (Dialog(tab="Storage"));
  parameter Boolean includeHeatTransfer=false "false for neglecting heat transfer" annotation (Dialog(tab="Storage"));
  parameter Modelica.Units.SI.Volume V_geo=1e5 "Geometric volume of storage" annotation (Dialog(tab="Storage"));
  parameter Modelica.Units.SI.Height height=3.779*V_geo^(1/3) "Height of storage" annotation (Dialog(tab="Storage"));
  parameter Modelica.Units.SI.CoefficientOfHeatTransfer alpha_nom=4 "Heat transfer coefficient inside the storage cylinder" annotation (Dialog(tab="Storage"));
  parameter Modelica.Units.SI.Mass m_start=1 "Stored gas mass at t=0" annotation (Dialog(tab="Storage"));
  parameter Modelica.Units.SI.Pressure p_start=simCenter.p_amb_const + simCenter.p_eff_2  "Pressure in storage at t=0" annotation (Dialog(tab="Storage")); // p_start = 17 bar
  parameter Modelica.Units.SI.ThermodynamicTemperature T_start=283.15 "Temperature of gas in storage at t=0" annotation (Dialog(tab="Storage"));
  parameter Modelica.Units.SI.Pressure p_out = 30e5;
  parameter Modelica.Units.SI.Pressure p_maxLow=p_maxHigh - 1e5 "Lower limit of the maximum pressure in storage" annotation (Dialog(tab="Storage", group="Control"));
  parameter Modelica.Units.SI.Pressure p_maxHigh=p_out "Upper limit of the maximum pressure in storage" annotation (Dialog(tab="Storage", group="Control"));
  parameter Modelica.Units.SI.Temperature T_out=283.15 "Hydrogen output temperature" annotation(Dialog(group="Fundamental Definitions"));
  parameter Modelica.Units.SI.Efficiency eta_n(
    min=0,
    max=1)=0.655 "Nominal efficency refering to the GCV (min = 0, max = 1)" annotation (Dialog(tab="General", group="Electrolyzer"));

  TransiEnt.Basics.Interfaces.Gas.RealGasPortIn H2PortIn(Medium=medium) "inlet flow" annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
  TransiEnt.Basics.Interfaces.Gas.RealGasPortOut H2PortOut(Medium=medium) "outlet flow" annotation (Placement(transformation(extent={{90,-10},{110,10}})));
  TransiEnt.Basics.Interfaces.General.PressureOut pressureTank annotation (Placement(transformation(extent={{96,44},{116,64}})));
  StorageSystem.TankSOC tankSOC annotation (Placement(transformation(extent={{22,66},{42,86}})));
  Modelica.Blocks.Interfaces.RealOutput socTank annotation (Placement(transformation(extent={{96,66},{116,86}})));

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
equation
  connect(storage.gasPortOut, H2PortOut) annotation (Line(
      points={{9.45,1},{9.45,0},{100,0}},
      color={255,255,0},
      thickness=1.5));
  connect(storage.p_gas, pressureTank) annotation (Line(points={{0,8},{0,54},{106,54}},    color={0,0,127}));
  connect(tankSOC.currentTankPressure, storage.p_gas) annotation (Line(points={{21.8,75.8},{0,75.8},{0,8}},    color={0,0,127}));
  connect(tankSOC.tankSOC, socTank) annotation (Line(points={{42.6,76},{106,76}}, color={0,0,127}));
  connect(H2PortIn, storage.gasPortIn) annotation (Line(
      points={{-100,0},{-20,0},{-20,1},{-7.35,1}},
      color={255,255,0},
      thickness=1.5));
  connect(H2PortOut, H2PortOut) annotation (Line(
      points={{100,0},{100,0}},
      color={255,255,0},
      thickness=1.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-54,82},{66,62}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid),
        Line(points={{-54,72},{-54,-68}}, color={28,108,200}),
        Ellipse(
          extent={{-54,-58},{66,-78}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid),
        Line(points={{66,72},{66,-68}}, color={28,108,200})}), Diagram(coordinateSystem(preserveAspectRatio=false)));
end H2StorageSystem_nonCompressed;
