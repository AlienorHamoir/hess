within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Tests;
model Test_PEMElectrolyzerL2_nonCompStorage "Test of PEM Electrolyzer L2 connection to non-compressed storage, cooling model and parameters calibration"

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
  Modelica.Blocks.Sources.Ramp PowerRamp(
    offset=0,
    startTime=0,
    duration=500,
    height=4e3)   annotation (Placement(transformation(extent={{-36,38},{-16,58}})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{14,74},{34,94}})));
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

  inner Modelica.Fluid.System system(energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics=Modelica.Fluid.Types.Dynamics.SteadyStateInitial)     annotation (
    Placement(visible = true, transformation(origin={84,84},      extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  StorageSystem.H2StorageSystem_nonCompressed H2StorageSystem(
    start_pressure=true,
    includeHeatTransfer=false,
    eta_n=0.75,
    V_geo=1,
    p_out=2070000,
    p_start=1700000) annotation (Placement(transformation(extent={{10,-70},{30,-50}})));
  Electrolyzer.Systems.SystemElectrolyzerL2_nonCompStorage ElectrolyzerSystem(
    electrolyzer(temperature(
        k_p=0.5,
        tau_i=5e-7,
        tau_d=0.1,
        N_d=0.6,
        N_i=10,
        PID_T_max(y=323.15))),
    usePowerPort=true,
    medium=medium,
    medium_coolant=medium_coolant,
    t_overload=900,
    m_flow_start=1e-4,
    P_el_min=275,
    p_out=2070000,
    T_out(displayUnit="K"),
    useHeatPort=true,
    useFluidCoolantPort=false,
    T_out_coolant_target=323.15,
    externalMassFlowControl=false) annotation (Placement(transformation(extent={{-14,-14},{14,14}})));
equation
  connect(MassflowRamp.y, H2massSink.m_flow) annotation (Line(points={{83,-28},{88,-28},{88,-66},{70,-66}}, color={0,0,127}));
  connect(H2StorageSystem.H2PortOut, H2massSink.gasPort) annotation (Line(
      points={{30,-60},{48,-60}},
      color={255,255,0},
      thickness=1.5));
  connect(ElectricGrid_0thOrder.epp, ElectrolyzerSystem.epp) annotation (Line(
      points={{-54,0},{-14,0}},
      color={0,135,135},
      thickness=0.5));
  connect(ElectrolyzerSystem.P_el_set, PowerRamp.y) annotation (Line(points={{0,14.56},{0,48},{-15,48}}, color={0,127,127}));
  connect(ElectrolyzerSystem.gasPortOut, H2StorageSystem.H2PortIn) annotation (Line(
      points={{0,-13.86},{0,-60},{10,-60}},
      color={255,255,0},
      thickness=1.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=10000,
      Interval=1,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"));
end Test_PEMElectrolyzerL2_nonCompStorage;
