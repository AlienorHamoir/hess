within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Tests;
model Test_PEMElectrolyzerL2_nonCompressedStorage "Test of PEM Electrolyzer L2 connection to non-compressed storage, cooling model and parameters calibration"

  extends TransiEnt.Basics.Icons.Checkmodel;
  import TransiEnt;
  import Modelica.Units.SI;

   // Import necessary packages for .csv reading
  import Modelica.Utilities.Streams.readRealMatrix;
  import Modelica.Utilities.Files.loadResource;


  // Medium declaration
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_coolant=simCenter.fluid1;

  TransiEnt.Components.Boundaries.Electrical.ActivePower.Frequency ElectricGrid_0thOrder annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-64,-12})));
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
    eta_n=0.851,
    V_geo=1.6,
    p_out=3000000,
    p_start=2000000) annotation (Placement(transformation(extent={{10,-70},{30,-50}})));
  Electrolyzer.Systems.SystemElectrolyzerL2_nonCompressedStorage ElectrolyzerSystem(
    electrolyzer(temperature(
        k_p=1000,
        tau_i=0.18,
        PID_T_max(y=328.15))),
    electrolyzer(voltage(humidity_const=21)),
    electrolyzer(massFlow(eta_F=1)),
    electrolyzer(pressure(p_mem_grad=17.1e5)),
    electrolyzer(T_op_start=50+273.15),
    usePowerPort=true,
    medium=medium,
    medium_coolant=medium_coolant,
    m_flow_start=1e-4,
    P_el_min=500,
    p_out=3000000,
    T_out(displayUnit="K"),
    useHeatPort=true,
    useFluidCoolantPort=false,
    T_out_coolant_target=323.15,
    externalMassFlowControl=false) annotation (Placement(transformation(extent={{-14,-26},{14,2}})));
    //     electrolyzer(T_op_start=50+273.15), // we use this line when we use Pstat as power setpoints
  Modelica.Blocks.Sources.Ramp PowerRampOperating(
    offset=0,
    startTime=0,
    duration=13000,
    height=4.725e3)
                  "Based on operating power curve- use P_el_set = P_el" annotation (Placement(transformation(extent={{-56,64},{-36,84}})));
  Modelica.Blocks.Sources.Ramp PowerRampCharacterization(
    offset=500,
    startTime=0,
    duration=13000,
    height=4.1e3) "Based on characterization power curve - use P_el_set = P_el" annotation (Placement(transformation(extent={{-56,26},{-36,46}})));
  Modelica.Blocks.Sources.Ramp PowerRampTest(
    offset=0,
    startTime=200,
    duration=1000,
    height=5e3) "Random power curve - use P_el_set = P_el_tot" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={52,52})));
  Modelica.Blocks.Sources.Ramp PowerRampCharacterization1(
    offset=3150,
    startTime=0,
    duration=4500,
    height=1450)  "Based on characterization power curve - use P_el_set = P_el" annotation (Placement(transformation(extent={{-90,26},{-70,46}})));
  Modelica.Blocks.Sources.Ramp PowerRampOperating1(
    offset=1475,
    startTime=100,
    duration=9000,
    height=3250)  "Based on operating power curve- use P_el_set = P_el" annotation (Placement(transformation(extent={{-90,64},{-70,84}})));


    // Path to the CSV file
  parameter String filePath = Modelica.Utilities.Files.loadResource("file:///C:/Users/alienor/Documents/hess/H2Microgrid_TransiEnt/Resources/Pdata.csv");

  Modelica.Blocks.Sources.CombiTimeTable Pdata(
    tableOnFile=true,
    tableName="Pdata",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/Pdata.txt"),
    verboseRead=true,
    columns={2,3,4},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=0.2) "Includes Pdata, Udata, Idata;  stop time is 4674 (timestep is 0.2 ms)"
                                                                     annotation (Placement(transformation(extent={{-26,26},{-6,46}})));
  Modelica.Blocks.Sources.CombiTimeTable Pstat(
    tableOnFile=true,
    tableName="Pstat",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/Pstat.txt"),
    verboseRead=true,
    columns={2,3,4},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=0.2) "Includes Pstat, Ustat, Istat; T_op_start must be 50 degC and stop time is 2280 sec (timestep is 0.2 ms)"                          annotation (Placement(transformation(extent={{-26,64},{-6,84}})));
  Modelica.Blocks.Sources.CombiTimeTable TempPressure(
    tableOnFile=true,
    tableName="TPdata",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/TPdata.txt"),
    verboseRead=true,
    columns={2,3},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=1) "Includes experimental H2 output Temperature and Stack Pressure over Pdata;  stop time is 23372" annotation (Placement(transformation(extent={{-88,-80},{-68,-60}})));
  Modelica.Blocks.Sources.CombiTimeTable StairSignal(table=[0,0; 499,0; 500,500; 999,500; 1000,1000; 1499,1000; 1500,1500; 1999,1500; 2000,2000; 2499,2000; 2500,2500; 2999,2500; 3000,3000; 3499,3000; 3500,3500; 3999,3500; 4000,4000; 4499,4000; 4500,4500; 4999,4500; 5000,5000; 5500,5000], tableOnFile=false) "create a stair-step signal for efficiency computation" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={52,22})));
equation
  connect(MassflowRamp.y, H2massSink.m_flow) annotation (Line(points={{83,-28},{88,-28},{88,-66},{70,-66}}, color={0,0,127}));
  connect(H2StorageSystem.H2PortOut, H2massSink.gasPort) annotation (Line(
      points={{30,-60},{48,-60}},
      color={255,255,0},
      thickness=1.5));
  connect(ElectricGrid_0thOrder.epp, ElectrolyzerSystem.epp) annotation (Line(
      points={{-54,-12},{-14,-12}},
      color={0,135,135},
      thickness=0.5));
  connect(ElectrolyzerSystem.gasPortOut, H2StorageSystem.H2PortIn) annotation (Line(
      points={{0,-25.86},{0,-60},{10,-60}},
      color={255,255,0},
      thickness=1.5));
  connect(PowerRampTest.y, ElectrolyzerSystem.P_el_set) annotation (Line(points={{41,52},{0,52},{0,2.56}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=4200,
      Interval=1,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"));
end Test_PEMElectrolyzerL2_nonCompressedStorage;
