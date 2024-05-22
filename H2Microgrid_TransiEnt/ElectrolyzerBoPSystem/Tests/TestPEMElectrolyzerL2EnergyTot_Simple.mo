within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Tests;
model TestPEMElectrolyzerL2EnergyTot_Simple "Test of PEM Electrolyzer L2 connection to storage"

  extends TransiEnt.Basics.Icons.Checkmodel;
  import TransiEnt;
  import Modelica.Units.SI;

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;
//*** DEFINE REPLACEABLE PACKAGES ***//
  // Medium declaration
  replaceable package medium_coolant = Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby Modelica.Media.Interfaces.PartialMedium;
  // Parameter definition
  parameter Real m_system_coolant(unit = "kg") = 44 "Coolant system mass";
    // Physical parameters
  parameter Real mass(unit = "kg") = 1 "Mass of the cell";
  parameter Real volume(unit = "m3") = 0.001 "Volume of the cell";
  // Thermal parameters
  parameter Real heatCapacity(unit = "J/(kg.K)") = 800 "Specific Heat Capacity";
//protected
//   function plotResult
//
//   constant String resultFileName = "TestPEMElectrolyzer_L1_Dynamics.mat";
//
//   output String resultFile;
//
//   algorithm
//     clearlog();
//     assert(cd(Modelica.Utilities.System.getEnvironmentVariable(TransiEnt.Basics.Types.WORKINGDIR)), "Error changing directory: Working directory must be set as environment variable with name 'workingdir' for this script to work.");
//     resultFile :=TransiEnt.Basics.Functions.fullPathName(Modelica.Utilities.System.getEnvironmentVariable(TransiEnt.Basics.Types.WORKINGDIR) + "/" + resultFileName);
//     removePlots(false);
//     createPlot(id=1, position={0, 0, 1563, 749}, y={"electrolyzer_0thOrder.dynamics.H_flow_H2", "electrolyzer_1stOrder.dynamics.H_flow_H2","electrolyzer_2ndOrder.dynamics.H_flow_H2"},range={0.0,52.0,-100000.0,1100000.0},erase=false,grid=true,filename=resultFile,colors={{28,108,200},{238,46,47},{0,140,72}});
//   end plotResult;

  TransiEnt.Components.Boundaries.Electrical.ActivePower.Frequency ElectricGrid_0thOrder annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-64,0})));
  Modelica.Blocks.Sources.Ramp PowerRamp(
    offset=0,
    startTime=0,
    duration=500,
    height=5e3)   annotation (Placement(transformation(extent={{-86,40},{-66,60}})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{-88,-96},{-68,-76}})));
  StorageSystem.H2StorageSystem_Compressed h2StorageSystem(
    start_pressure=true,
    includeHeatTransfer=false,
    eta_n=0.75,
    V_geo=1,
    p_out=20e5,
    p_start=17e5)  annotation (Placement(transformation(extent={{12,-70},{32,-50}})));
  inner TransiEnt.ModelStatistics modelStatistics annotation (Placement(transformation(extent={{-54,-96},{-34,-76}})));

//      integrateMassFlow=true

  Modelica.Blocks.Sources.Ramp MassflowRamp(
    offset=0,
    startTime=20,
    duration=60,
    height=1e-5) annotation (Placement(transformation(extent={{62,-38},{82,-18}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_Txim_flow boundary_Txim_flow(medium=medium, variable_m_flow=true) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={58,-62})));
  Electrolyzer.SystemElectrolyzerL2_Simple_EnergyTot systemElectrolyzerL2_Simple_EnergyTot(
    usePowerPort=true,
    medium=medium,
    t_overload=900,
    m_flow_start=1e-4,
    P_el_min=275,
    k=1e11,
    p_out=20e5,
    useHeatPort=true,
    electrolyzer(
      useFluidCoolantPort=true,
      useHeatPort=false,
      externalMassFlowControl=true))
                      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

  inner Modelica.Fluid.System system(energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics=Modelica.Fluid.Types.Dynamics.SteadyStateInitial)     annotation (
    Placement(visible = true, transformation(origin={-16,-86},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(h2StorageSystem.H2PortOut, boundary_Txim_flow.gasPort) annotation (Line(
      points={{32.4,-62.6},{44,-62.6},{44,-62},{48,-62}},
      color={255,255,0},
      thickness=1.5));
  connect(MassflowRamp.y, boundary_Txim_flow.m_flow) annotation (Line(points={{83,-28},{88,-28},{88,-68},{70,-68}}, color={0,0,127}));
  connect(systemElectrolyzerL2_Simple_EnergyTot.epp, ElectricGrid_0thOrder.epp) annotation (Line(
      points={{-10,0},{-54,0}},
      color={0,135,135},
      thickness=0.5));
  connect(systemElectrolyzerL2_Simple_EnergyTot.gasPortOut, h2StorageSystem.H2PortIn) annotation (Line(
      points={{0,-9.9},{0,-62.6},{12,-62.6}},
      color={255,255,0},
      thickness=1.5));
  connect(systemElectrolyzerL2_Simple_EnergyTot.P_el_set, PowerRamp.y) annotation (Line(points={{0,10.4},{0,50},{-65,50}}, color={0,127,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(Interval=1, __Dymola_Algorithm="Dassl"));
end TestPEMElectrolyzerL2EnergyTot_Simple;
