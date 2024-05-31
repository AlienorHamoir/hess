within H2Microgrid_TransiEnt;
model HESS
  ElectrolyzerBoPSystem.Electrolyzer.Systems.SystemElectrolyzerL2_compStorage systemElectrolyzerL2_Simple_EnergyTot_Inv_Cooling annotation (Placement(transformation(extent={{-60,-80},{-20,-40}})));
  StorageSystem.H2StorageSystem_Compressed h2StorageSystem_Compressed annotation (Placement(transformation(extent={{20,-40},{60,20}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn P_set_electrolyzer "Electric power input" annotation (Placement(transformation(extent={{-112,-42},{-92,-22}})));
  Modelica.Blocks.Interfaces.RealOutput socTankH2 annotation (Placement(transformation(extent={{96,-10},{116,10}})));
equation
  connect(systemElectrolyzerL2_Simple_EnergyTot_Inv_Cooling.gasPortOut, h2StorageSystem_Compressed.H2PortIn) annotation (Line(
      points={{-40,-79.8},{-40,-88},{40.2,-88},{40.2,-40.3}},
      color={255,255,0},
      thickness=1.5));
  connect(systemElectrolyzerL2_Simple_EnergyTot_Inv_Cooling.P_el_set, P_set_electrolyzer) annotation (Line(points={{-40,-39.2},{-40,-32},{-102,-32}}, color={0,127,127}));
  connect(h2StorageSystem_Compressed.socTank, socTankH2) annotation (Line(points={{61.2,15.8},{92,15.8},{92,0},{106,0}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={255,213,170},
          fillPattern=FillPattern.Solid), Text(
          extent={{-90,42},{90,-36}},
          textColor={0,0,0},
          textString="HESS")}),                                  Diagram(coordinateSystem(preserveAspectRatio=false)));
end HESS;
