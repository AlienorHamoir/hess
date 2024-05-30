within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Tests;
model TestPEMElectrolyzerL2EnergyTot_Inv_Cooling "Test of PEM Electrolyzer L2 connection to storage"

  extends TransiEnt.Basics.Icons.Checkmodel;
  import TransiEnt;
  import Modelica.Units.SI;

  // Medium declaration
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_coolant=simCenter.fluid1;
  // Parameter definition
//   parameter Real m_system_coolant(unit = "kg") = 44 "Coolant system mass";
//     // Physical parameters
//   parameter Real mass(unit = "kg") = 1 "Mass of the cell";
//   parameter Real volume(unit = "m3") = 0.001 "Volume of the cell";
//   // Thermal parameters
//   parameter Real heatCapacity(unit = "J/(kg.K)") = 800 "Specific Heat Capacity";

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
  StorageSystem.H2StorageSystem_Compressed h2StorageSystem(
    start_pressure=true,
    includeHeatTransfer=false,
    eta_n=0.75,
    V_geo=1,
    p_out=2000000,
    p_start=1700000)
                   annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={24,-60})));
  inner TransiEnt.ModelStatistics modelStatistics annotation (Placement(transformation(extent={{46,74},{66,94}})));

  Modelica.Blocks.Sources.Ramp MassflowRamp(
    offset=0,
    startTime=20,
    duration=60,
    height=1e-5) annotation (Placement(transformation(extent={{62,-38},{82,-18}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_Txim_flow boundary_Txim_flow(medium=medium, variable_m_flow=true) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={58,-60})));

  inner Modelica.Fluid.System system(energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics=Modelica.Fluid.Types.Dynamics.SteadyStateInitial)     annotation (
    Placement(visible = true, transformation(origin={84,84},      extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Electrolyzer.SystemElectrolyzerL2_Simple_EnergyTot_Inv_Cooling systemElectrolyzerL2_Simple_EnergyTot_Inv_Cooling(
   usePowerPort=true,
    medium=medium,
    medium_coolant=medium_coolant,
    t_overload=900,
    m_flow_start=1e-4,
    P_el_min=275,
    p_out=2000000,
    useHeatPort=true,
    useFluidCoolantPort=false,
    T_out_coolant_target=323.15,
    externalMassFlowControl=false,
    electrolyzer(temperature(
        k_p=0.5, tau_i=5e-7,
        tau_d=0.1, N_d=0.6,  N_i=10,           PID_T_max(y=323.15))))
                                   annotation (Placement(transformation(extent={{-14,-14},{14,14}})));
equation
  connect(h2StorageSystem.H2PortOut, boundary_Txim_flow.gasPort) annotation (Line(
      points={{33.8,-60},{48,-60}},
      color={255,255,0},
      thickness=1.5));
  connect(MassflowRamp.y, boundary_Txim_flow.m_flow) annotation (Line(points={{83,-28},{88,-28},{88,-66},{70,-66}}, color={0,0,127}));
  connect(systemElectrolyzerL2_Simple_EnergyTot_Inv_Cooling.epp, ElectricGrid_0thOrder.epp) annotation (Line(
      points={{-14,0},{-54,0}},
      color={0,135,135},
      thickness=0.5));
  connect(systemElectrolyzerL2_Simple_EnergyTot_Inv_Cooling.P_el_set, PowerRamp.y) annotation (Line(points={{0,14.56},{0,48},{-15,48}},                       color={0,127,127}));
  connect(systemElectrolyzerL2_Simple_EnergyTot_Inv_Cooling.gasPortOut, h2StorageSystem.H2PortIn) annotation (Line(
      points={{0,-13.86},{0,-60.1},{13.9,-60.1}},
      color={255,255,0},
      thickness=1.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=10000,
      Interval=1,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"));
end TestPEMElectrolyzerL2EnergyTot_Inv_Cooling;
