within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Tests;
model TestStorage_DryingSystem

  extends TransiEnt.Basics.Icons.Checkmodel;
  import TransiEnt;
  import      Modelica.Units.SI;

public
  parameter SI.Power P_el_n=1e6 "Nominal electrical power of the electrolyzer";
  parameter SI.Power P_el_min=0.05*P_el_n "Minimal electrical power of the electrolyzer";
  parameter SI.Power P_el_max=1.0*P_el_n "Maximal electrical power of the electrolyzer";
  parameter SI.Temperature T_out=273.15+15 "Temperature of the produced hydrogen";
  parameter SI.Efficiency eta_n=0.75 "Nominal efficiency of the electrolyzer";
  parameter SI.Pressure p_out=50e5 "Pressure of the produced hydrogen";

  parameter TILMedia.VLEFluidTypes.TILMedia_SplineWater water;
  parameter TransiEnt.Basics.Media.Gases.VLE_VDIWA_SG4_var medium_gas;


  inner TransiEnt.SimCenter
                  simCenter annotation (Placement(transformation(extent={{-86,74},{-66,94}})));
  Modelica.Blocks.Sources.Ramp ramp2(
    offset=30e5,
    duration=1000,
    startTime=7000,
    height=-20e5)  annotation (Placement(transformation(extent={{90,64},{70,84}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_Txim_flow source(
    medium=medium_gas,
    variable_m_flow=false,
    variable_T=false,
    variable_xi=false)
                      annotation (Placement(transformation(extent={{-68,10},{-48,30}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_pTxi sink_syngas(
    medium=medium_gas,
    variable_p=false,
    p_const=100000) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={52,20})));
  StorageSystem.H2StorageSystem_Compressed h2StorageSystem annotation (Placement(transformation(extent={{0,12},{20,32}})));
  inner TransiEnt.ModelStatistics modelStatistics annotation (Placement(transformation(extent={{-96,-82},{-76,-62}})));
  H2DryingSystem.H2DryingSubsystem h2DryingSubsystem(water=water, medium_gas=medium_gas) annotation (Placement(transformation(extent={{-34,12},{-14,32}})));
equation
  connect(h2StorageSystem.H2PortOut, sink_syngas.gasPort) annotation (Line(
      points={{20.4,20.6},{36,20.6},{36,20},{42,20}},
      color={255,255,0},
      thickness=1.5));
  connect(source.gasPort, h2DryingSubsystem.gasPortIn) annotation (Line(
      points={{-48,20},{-48,22},{-33.8,22}},
      color={255,255,0},
      thickness=1.5));
  connect(h2DryingSubsystem.gasPortOut, h2StorageSystem.H2PortIn) annotation (Line(
      points={{-14,22},{-14,20.6},{0,20.6}},
      color={255,255,0},
      thickness=1.5));
  connect(source.m_flow, ramp2.y) annotation (Line(points={{-70,26},{-70,62},{64,62},{64,74},{69,74}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end TestStorage_DryingSystem;
