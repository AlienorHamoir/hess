within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Tests;
model TestPEMElectrolyzerL1_CoolingSystem

  extends TransiEnt.Basics.Icons.Checkmodel;

  ElectrolyzerBoPSystem.Electrolyzer.PEMElectrolyzer_L1 electrolyzer(
    useFluidCoolantPort=true,
    externalMassFlowControl=true,
    T_out_coolant_target=333.15,
    heatFlow_externalMassFlowControl(Medium=simCenter.fluid1)) annotation (Placement(transformation(extent={{-52,20},{-32,40}})));
  TransiEnt.Components.Boundaries.Electrical.ActivePower.Frequency ElectricGrid_0thOrder annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-78,30})));
  Modelica.Blocks.Sources.Ramp PowerRamp(
    offset=2,
    startTime=0,
    duration=3,
    height=100000) annotation (Placement(transformation(extent={{-78,70},{-58,90}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_pTxi boundary_pTxi(medium=simCenter.gasModel3) annotation (Placement(transformation(extent={{6,20},{-14,40}})));
  Buildings.Experimental.DHC.EnergyTransferStations.BaseClasses.Pump_m_flow pump_m_flow annotation (Placement(transformation(extent={{26,-26},{46,-6}})));
  Buildings.Experimental.DHC.Loads.HotWater.BaseClasses.HeatExchangerPumpController conPum annotation (Placement(transformation(extent={{-36,-82},{-16,-62}})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{-88,-96},{-68,-76}})));
  TransiEnt.Components.Turbogroups.Hydroturbine hydroturbine annotation (Placement(transformation(extent={{-40,-28},{-20,-8}})));
equation
  connect(ElectricGrid_0thOrder.epp, electrolyzer.epp) annotation (Line(
      points={{-68,30},{-52,30}},
      color={0,135,135},
      thickness=0.5));
  connect(PowerRamp.y, electrolyzer.P_el_set) annotation (Line(points={{-57,80},{-46,80},{-46,42}}, color={0,0,127}));
  connect(boundary_pTxi.gasPort, electrolyzer.gasPortOut) annotation (Line(
      points={{-14,30},{-32,30}},
      color={255,255,0},
      thickness=1.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end TestPEMElectrolyzerL1_CoolingSystem;
