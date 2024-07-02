within H2Microgrid_TransiEnt;
model HESS_nonCompressed "HESS without compressed storage"

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;

  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn P_set_electrolyzer "Electric power input" annotation (Placement(transformation(extent={{-120,-70},{-100,-50}})));
  Modelica.Blocks.Interfaces.RealOutput socTankH2 annotation (Placement(transformation(extent={{96,-10},{116,10}})));
  FuelCellBoPSystem.FuelCell.SystemPEMFC systemPEMFC annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-40,60})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn P_set_FC "Input for FC power production setpoint" annotation (Placement(transformation(extent={{-114,54},{-94,74}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_FC "Actual FC electrical power production" annotation (Placement(transformation(extent={{96,72},{116,92}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_electrolyzer annotation (Placement(transformation(extent={{98,-78},{118,-58}})));
  StorageSystem.H2StorageSystem_nonCompressed h2StorageSystem_nonCompressed annotation (Placement(transformation(extent={{20,-50},{56,6}})));
  ElectrolyzerBoPSystem.Electrolyzer.Systems.SystemElectrolyzerL2_nonCompressedStorage systemElectrolyzerL2_nonCompressedStorage(usePowerPort=false)
                                                                                                                                 annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-40,-60})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{-64,0},{-44,20}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_Txim_flow H2massSink(medium=medium, variable_m_flow=true) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={40,44})));
equation
  connect(systemPEMFC.P_el_set, P_set_FC) annotation (Line(points={{-61.6,64.8},{-90,64.8},{-90,64},{-104,64}}, color={0,127,127}));
  connect(systemPEMFC.P_FC_tot, P_FC) annotation (Line(
      points={{-32,80.8},{106,80.8},{106,82}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(systemElectrolyzerL2_nonCompressedStorage.P_el_set, P_set_electrolyzer) annotation (Line(points={{-60.8,-60},{-110,-60}}, color={0,127,127}));
  connect(systemElectrolyzerL2_nonCompressedStorage.electrolyzerPowerOut, P_electrolyzer) annotation (Line(
      points={{-24.8,-39.2},{-24.8,-32},{14,-32},{14,-68},{108,-68}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(systemElectrolyzerL2_nonCompressedStorage.gasPortOut, h2StorageSystem_nonCompressed.H2PortIn) annotation (Line(
      points={{-20.2,-60},{16,-60},{16,-22},{20,-22}},
      color={255,255,0},
      thickness=1.5));
  connect(h2StorageSystem_nonCompressed.socTank, socTankH2) annotation (Line(points={{57.08,-0.72},{92,-0.72},{92,0},{106,0}}, color={0,0,127}));
  connect(systemPEMFC.mflowH2_FC, H2massSink.m_flow) annotation (Line(points={{-18,48.8},{20,48.8},{20,50},{28,50}}, color={0,0,127}));
  connect(H2massSink.gasPort, h2StorageSystem_nonCompressed.H2PortOut) annotation (Line(
      points={{50,44},{62,44},{62,-22},{56,-22}},
      color={255,255,0},
      thickness=1.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={255,213,170},
          fillPattern=FillPattern.Solid), Text(
          extent={{-90,42},{90,-36}},
          textColor={0,0,0},
          textString="HESS
30 bar")}),                                                      Diagram(coordinateSystem(preserveAspectRatio=false)));
end HESS_nonCompressed;
