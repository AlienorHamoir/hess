within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.OLD;
model TestPEMElectrolyzerL2EnergyTot_Inv_innerCooling "Test of PEM Electrolyzer L2 connection to storage"

  extends TransiEnt.Basics.Icons.Checkmodel;
  import TransiEnt;
  import Modelica.Units.SI;

  // Medium declaration
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_coolant=simCenter.fluid1;
  // Parameter definition
  parameter Real m_system_coolant(unit = "kg") = 44 "Coolant system mass";
    // Physical parameters
  parameter Real mass(unit = "kg") = 1 "Mass of the cell";
  parameter Real volume(unit = "m3") = 0.001 "Volume of the cell";
  // Thermal parameters
  parameter Real heatCapacity(unit = "J/(kg.K)") = 800 "Specific Heat Capacity";

  TransiEnt.Components.Boundaries.Electrical.ActivePower.Frequency ElectricGrid_0thOrder annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-64,0})));
  Modelica.Blocks.Sources.Ramp PowerRamp(
    offset=0,
    startTime=0,
    duration=500,
    height=5e3)   annotation (Placement(transformation(extent={{-36,38},{-16,58}})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{14,74},{34,94}})));
  StorageSystem.H2StorageSystem_Compressed h2StorageSystem(
    start_pressure=true,
    includeHeatTransfer=false,
    eta_n=0.75,
    V_geo=1,
    p_out=20e5,
    p_start=17e5)  annotation (Placement(transformation(extent={{14,-70},{34,-50}})));
  inner TransiEnt.ModelStatistics modelStatistics annotation (Placement(transformation(extent={{46,74},{66,94}})));

  Modelica.Blocks.Sources.Ramp MassflowRamp(
    offset=0,
    startTime=20,
    duration=60,
    height=1e-5) annotation (Placement(transformation(extent={{62,-38},{82,-18}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_Txim_flow boundary_Txim_flow(medium=medium, variable_m_flow=true) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={58,-62})));

  inner Modelica.Fluid.System system(energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics=Modelica.Fluid.Types.Dynamics.SteadyStateInitial)     annotation (
    Placement(visible = true, transformation(origin={84,84},      extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SystemElectrolyzerL2_Simple_EnergyTot_Inv_innerCooling systemElectrolyzerL2_Simple_EnergyTot_Inv(
    usePowerPort=true,
    medium=medium,
    medium_coolant=medium_coolant,
    t_overload=900,
    m_flow_start=1e-4,
    P_el_min=275,
    k=1e11,
    p_out=2000000,
    useHeatPort=false,
    useFluidCoolantPort=false,
    T_out_coolant_target=323.15,
    externalMassFlowControl=false,
    electrolyzer1(temperature(cooling_PID(
          controllerType=Modelica.Blocks.Types.SimpleController.PID,
          k=1,
          Tau_i=0.00001,
          Tau_d=0.01,
          Ni=0.9)))) annotation (Placement(transformation(extent={{-12,-16},{16,14}})));
equation
  connect(h2StorageSystem.H2PortOut, boundary_Txim_flow.gasPort) annotation (Line(
      points={{24,-50.2},{44,-50.2},{44,-62},{48,-62}},
      color={255,255,0},
      thickness=1.5));
  connect(MassflowRamp.y, boundary_Txim_flow.m_flow) annotation (Line(points={{83,-28},{88,-28},{88,-68},{70,-68}}, color={0,0,127}));
  connect(ElectricGrid_0thOrder.epp, systemElectrolyzerL2_Simple_EnergyTot_Inv.epp) annotation (Line(
      points={{-54,0},{-22,0},{-22,-1},{-12,-1}},
      color={0,135,135},
      thickness=0.5));
  connect(systemElectrolyzerL2_Simple_EnergyTot_Inv.P_el_set, PowerRamp.y) annotation (Line(points={{2,14.6},{2,48},{-15,48}},      color={0,127,127}));
  connect(systemElectrolyzerL2_Simple_EnergyTot_Inv.gasPortOut, h2StorageSystem.H2PortIn) annotation (Line(
      points={{2,-15.85},{2,-70.1},{24.1,-70.1}},
      color={255,255,0},
      thickness=1.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=10000,
      Interval=1,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"));
end TestPEMElectrolyzerL2EnergyTot_Inv_innerCooling;
