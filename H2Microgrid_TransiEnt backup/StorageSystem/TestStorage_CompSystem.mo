within H2Microgrid_TransiEnt.StorageSystem;
model TestStorage_CompSystem

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
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;

  inner TransiEnt.SimCenter
                  simCenter annotation (Placement(transformation(extent={{-86,74},{-66,94}})));
  Modelica.Blocks.Sources.Ramp ramp2(
    offset=30e5,
    duration=1000,
    startTime=7000,
    height=-20e5)  annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-74,50})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_Txim_flow source(
    medium=medium,
    variable_m_flow=false,
    variable_T=false,
    variable_xi=false)
                      annotation (Placement(transformation(extent={{-46,10},{-26,30}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_pTxi sink_syngas(
    medium=medium, p_const=1000000)
                  annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={52,20})));
  inner TransiEnt.ModelStatistics modelStatistics annotation (Placement(transformation(extent={{-96,-82},{-76,-62}})));
  H2StorageSystem_Compressed h2StorageSystem_Compressed(medium=medium) annotation (Placement(transformation(extent={{-4,12},{16,32}})));
equation
  connect(source.m_flow, ramp2.y) annotation (Line(points={{-48,26},{-50,26},{-50,50},{-63,50}},       color={0,0,127}));
  connect(source.gasPort, h2StorageSystem_Compressed.H2PortIn) annotation (Line(
      points={{-26,20},{-6,20},{-6,19.4},{-4,19.4}},
      color={255,255,0},
      thickness=1.5));
  connect(h2StorageSystem_Compressed.H2PortOut, sink_syngas.gasPort) annotation (Line(
      points={{16.4,19.4},{36,19.4},{36,20},{42,20}},
      color={255,255,0},
      thickness=1.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end TestStorage_CompSystem;
