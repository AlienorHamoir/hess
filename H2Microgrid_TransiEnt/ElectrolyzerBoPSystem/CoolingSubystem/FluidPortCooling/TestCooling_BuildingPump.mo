within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.CoolingSubystem.FluidPortCooling;
model TestCooling_BuildingPump
  extends TransiEnt.Basics.Icons.Checkmodel;
  import TransiEnt;
  import Modelica.Units.SI;

   parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_coolant=simCenter.fluid1;

  CoolingSubystem.FluidPortCooling.Cooling_BuildingPump cooling_BuildingPump annotation (Placement(transformation(extent={{-50,-12},{-30,8}})));
  TransiEnt.Components.Boundaries.FluidFlow.BoundaryVLE_hxim_flow boundaryVLE_hxim_flow(medium=medium_coolant) annotation (Placement(transformation(extent={{-2,10},{18,30}})));
  TransiEnt.Components.Boundaries.FluidFlow.BoundaryVLE_phxi boundaryVLE_phxi(medium=medium_coolant) annotation (Placement(transformation(extent={{-74,-40},{-54,-20}})));
  Modelica.Blocks.Sources.Ramp ramp(
    height=50,
    duration=100,
    startTime=5) annotation (Placement(transformation(extent={{-64,42},{-44,62}})));
  inner SimCenter simCenter annotation (Placement(transformation(extent={{54,-88},{74,-68}})));
equation
  connect(cooling_BuildingPump.fluidPortOut, boundaryVLE_phxi.fluidPortIn) annotation (Line(
      points={{-30,-10},{-26,-10},{-26,-30},{-54,-30}},
      color={175,0,0},
      thickness=0.5));
  connect(cooling_BuildingPump.fluidPortIn, boundaryVLE_hxim_flow.fluidPortOut) annotation (Line(
      points={{-30,-6},{24,-6},{24,20},{18,20}},
      color={175,0,0},
      thickness=0.5));
  connect(ramp.y, boundaryVLE_hxim_flow.m_flow) annotation (Line(points={{-43,52},{-12,52},{-12,26},{-4,26}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end TestCooling_BuildingPump;
