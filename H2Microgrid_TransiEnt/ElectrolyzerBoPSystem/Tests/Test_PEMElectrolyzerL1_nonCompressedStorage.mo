within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Tests;
model Test_PEMElectrolyzerL1_nonCompressedStorage "Test of PEM Electrolyzer L1 connection to non-compressed storage and parameters calibration"

  extends TransiEnt.Basics.Icons.Checkmodel;
  import TransiEnt;
  import Modelica.Units.SI;

  // Medium declaration
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_coolant=simCenter.fluid1;

  TransiEnt.Components.Boundaries.Electrical.ActivePower.Frequency ElectricGrid_0thOrder annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-64,0})));
  Modelica.Blocks.Sources.Ramp PowerRampCharacterization(
    offset=3000,
    startTime=0,
    duration=4500,
    height=1.6e3) "Based on characterization power curve - use P_el_set = P_el" annotation (Placement(transformation(extent={{-36,34},{-16,54}})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{14,74},{34,94}})));
  StorageSystem.H2StorageSystem_Compressed H2StorageSystem(
    start_pressure=true,
    includeHeatTransfer=false,
    eta_n=0.851,
    V_geo=1.6,
    p_out=3000000,
    p_start=2000000)
                   annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={24,-60})));
  inner TransiEnt.ModelStatistics modelStatistics annotation (Placement(transformation(extent={{46,74},{66,94}})));

  Modelica.Blocks.Sources.Ramp MassflowRamp(
    offset=0,
    startTime=20,
    duration=60,
    height=1e-5) annotation (Placement(transformation(extent={{62,-38},{82,-18}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_Txim_flow H2massSink(medium=medium, variable_m_flow=true) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={58,-60})));

  Modelica.Blocks.Sources.Ramp PowerRampOperating(
    offset=1500,
    startTime=0,
    duration=9000,
    height=3.2e3) "Based on operating power curve- use P_el_set = P_el" annotation (Placement(transformation(extent={{-36,68},{-16,88}})));
  Modelica.Blocks.Sources.Ramp PowerRampTest(
    offset=0,
    startTime=0,
    duration=1000,
    height=5e3) "Random power curve - use P_el_set = P_el_tot" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={40,40})));
  Electrolyzer.Systems.SystemElectrolyzerL1_nonCompressedStorage systemElectrolyzerL1_nonCompStorage(
    usePowerPort=true,
    medium=medium,
    medium_coolant=medium_coolant,
    m_flow_start=1e-5,
    P_el_min=500,
    p_out=3000000,
    useHeatPort=true,
    useFluidCoolantPort=false,
    T_out_coolant_target=323.15,
    externalMassFlowControl=false) annotation (Placement(transformation(extent={{-16,-16},{16,16}})));
equation
  connect(H2StorageSystem.H2PortOut, H2massSink.gasPort) annotation (Line(
      points={{33.8,-60},{48,-60}},
      color={255,255,0},
      thickness=1.5));
  connect(MassflowRamp.y, H2massSink.m_flow) annotation (Line(points={{83,-28},{88,-28},{88,-66},{70,-66}}, color={0,0,127}));
  connect(systemElectrolyzerL1_nonCompStorage.epp, ElectricGrid_0thOrder.epp) annotation (Line(
      points={{-16,0},{-54,0}},
      color={0,135,135},
      thickness=0.5));
  connect(systemElectrolyzerL1_nonCompStorage.P_el_set, PowerRampTest.y) annotation (Line(points={{0,16.64},{0,40},{29,40}},                       color={0,127,127}));
  connect(systemElectrolyzerL1_nonCompStorage.gasPortOut, H2StorageSystem.H2PortIn) annotation (Line(
      points={{0,-15.84},{0,-60.1},{13.9,-60.1}},
      color={255,255,0},
      thickness=1.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=3600,
      Interval=1,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"));
end Test_PEMElectrolyzerL1_nonCompressedStorage;
