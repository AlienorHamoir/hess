within ;
model PowerConnector_Example
  TransiEnt.Basics.Interfaces.Electrical.ActivePowerPort epp annotation (Placement(transformation(extent={{-80,-98},{-60,-78}})));
  TransiEnt.Components.Boundaries.Electrical.ActivePower.Frequency ElectricGrid_0thOrder annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-74,2})));
  TransiEnt.Consumer.Systems.PVBatteryPoolControl.PVBatteryPool pVBatteryPool(nLoadProfiles=1, nPVProfiles=1) annotation (Placement(transformation(extent={{14,-8},{34,12}})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{66,-88},{86,-68}})));
  Modelica.Blocks.Sources.Sine PVGen(
    amplitude=500,
    f=50,
    offset=100,
    startTime=10) annotation (Placement(transformation(extent={{-44,30},{-24,50}})));
  Modelica.Blocks.Sources.Trapezoid Load(
    amplitude=600,
    rising=40,
    width=100,
    falling=55,
    period=250,
    nperiod=-1,
    offset=30,
    startTime=10) annotation (Placement(transformation(extent={{-44,64},{-24,84}})));
equation
  connect(ElectricGrid_0thOrder.epp, pVBatteryPool.epp) annotation (Line(
      points={{-64,2},{4,2},{4,2.05},{14.1,2.05}},
      color={0,135,135},
      thickness=0.5));
  connect(PVGen.y, pVBatteryPool.P_el_PV[1]) annotation (Line(points={{-23,40},{4,40},{4,5.8},{13.6,5.8}}, color={0,0,127}));
  connect(Load.y, pVBatteryPool.P_el_load[1]) annotation (Line(points={{-23,74},{6,74},{6,9.8},{13.6,9.8}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    uses(TransiEnt(version="2.0.3"), Modelica(version="4.0.0")));
end PowerConnector_Example;
