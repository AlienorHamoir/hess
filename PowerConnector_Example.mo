within ;
model PowerConnector_Example
  TransiEnt.Basics.Interfaces.Electrical.ActivePowerPort epp annotation (Placement(transformation(extent={{-84,-62},{-64,-42}})));
  TransiEnt.Components.Boundaries.Electrical.ActivePower.Frequency ElectricGrid_0thOrder annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-74,2})));
  TransiEnt.Consumer.Systems.PVBatteryPoolControl.PVBatteryPool pVBatteryPool annotation (Placement(transformation(extent={{14,-8},{34,12}})));
equation
  connect(ElectricGrid_0thOrder.epp, pVBatteryPool.epp) annotation (Line(
      points={{-64,2},{4,2},{4,2.05},{14.1,2.05}},
      color={0,135,135},
      thickness=0.5));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    uses(TransiEnt(version="2.0.3")));
end PowerConnector_Example;
