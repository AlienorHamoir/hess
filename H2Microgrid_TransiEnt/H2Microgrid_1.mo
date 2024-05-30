within H2Microgrid_TransiEnt;
model H2Microgrid_1
  HESS hESS annotation (Placement(transformation(extent={{20,-60},{60,-20}})));
  SCooDER.Systems.SmartBuilding.smartBuilding_external smartBuilding_external(n_trans=1) annotation (Placement(transformation(extent={{-42,0},{-2,40}})));
  Modelica.Blocks.Sources.Constant Load(k=0) annotation (Placement(transformation(extent={{-92,4},{-72,24}})));
  Modelica.Blocks.Interfaces.RealInput P_charger1[20] annotation (Placement(transformation(extent={{-10,-74},{30,-34}})));
equation
  connect(Load.y, smartBuilding_external.P_building) annotation (Line(points={{-71,14},{-44,14}}, color={0,0,127}));
  connect(smartBuilding_external.P_charger, P_charger1) annotation (Line(points={{-44,8},{-52,8},{-52,-54},{10,-54}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end H2Microgrid_1;
