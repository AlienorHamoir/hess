within H2Microgrid_TransiEnt;
model HESS_Compressed
  ElectrolyzerBoPSystem.Electrolyzer.Systems.SystemElectrolyzerL2_CompressedStorage systemElectrolyzerL2_Simple_EnergyTot_Inv_Cooling annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-40,-60})));
  StorageSystem.H2StorageSystem_Compressed h2StorageSystem_Compressed annotation (Placement(transformation(extent={{20,-40},{60,20}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn P_set_electrolyzer "Electric power input" annotation (Placement(transformation(extent={{-120,-70},{-100,-50}})));
  Modelica.Blocks.Interfaces.RealOutput socTankH2 annotation (Placement(transformation(extent={{96,-10},{116,10}})));
  FuelCellBoPSystem.FuelCell.SystemPEMFC systemPEMFC annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-40,60})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn P_set_FC "Input for FC power production setpoint" annotation (Placement(transformation(extent={{-114,54},{-94,74}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_Txim_flow H2massSink(medium=medium, variable_m_flow=true) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={24,46})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_FC "Actual FC electrical power production" annotation (Placement(transformation(extent={{96,72},{116,92}})));
  TransiEnt.Components.Boundaries.Electrical.ActivePower.Frequency ElectricGrid_0thOrder annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-84,0})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_electrolyzer annotation (Placement(transformation(extent={{98,-78},{118,-58}})));
equation
  connect(systemElectrolyzerL2_Simple_EnergyTot_Inv_Cooling.gasPortOut, h2StorageSystem_Compressed.H2PortIn) annotation (Line(
      points={{-20.2,-60},{40.2,-60},{40.2,-40.3}},
      color={255,255,0},
      thickness=1.5));
  connect(systemElectrolyzerL2_Simple_EnergyTot_Inv_Cooling.P_el_set, P_set_electrolyzer) annotation (Line(points={{-60.8,-60},{-110,-60}}, color={0,127,127}));
  connect(h2StorageSystem_Compressed.socTank, socTankH2) annotation (Line(points={{61.2,13.4},{92,13.4},{92,0},{106,0}}, color={0,0,127}));
  connect(systemPEMFC.P_el_set, P_set_FC) annotation (Line(points={{-61.6,64.8},{-90,64.8},{-90,64},{-104,64}}, color={0,127,127}));
  connect(h2StorageSystem_Compressed.H2PortOut, H2massSink.gasPort) annotation (Line(
      points={{40,19.4},{40,46},{34,46}},
      color={255,255,0},
      thickness=1.5));
  connect(systemPEMFC.mflowH2_FC, H2massSink.m_flow) annotation (Line(points={{-18.2,51.8},{-18.2,52},{12,52}}, color={0,0,127}));
  connect(h2StorageSystem_Compressed.P_comp, systemElectrolyzerL2_Simple_EnergyTot_Inv_Cooling.CompressorPower) annotation (Line(
      points={{60.8,-20.8},{60.8,-26},{68,-26},{68,-90},{-28.4,-90},{-28.4,-81}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(systemPEMFC.P_FC_tot, P_FC) annotation (Line(
      points={{-32.8,82},{106,82},{106,82}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(systemPEMFC.epp, ElectricGrid_0thOrder.epp) annotation (Line(
      points={{-50.4,40},{-50,40},{-50,0},{-74,0}},
      color={0,135,135},
      thickness=0.5));
  connect(systemElectrolyzerL2_Simple_EnergyTot_Inv_Cooling.epp, ElectricGrid_0thOrder.epp) annotation (Line(
      points={{-40,-80},{-40,-88},{-68,-88},{-68,0},{-74,0}},
      color={0,135,135},
      thickness=0.5));
  connect(systemElectrolyzerL2_Simple_EnergyTot_Inv_Cooling.P_electrolyzer_tot, P_electrolyzer) annotation (Line(
      points={{-26,-38.8},{-8,-38.8},{-8,-68},{108,-68}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={255,213,170},
          fillPattern=FillPattern.Solid), Text(
          extent={{-90,42},{90,-36}},
          textColor={0,0,0},
          textString="HESS
350 bar")}),                                                     Diagram(coordinateSystem(preserveAspectRatio=false)));
end HESS_Compressed;
