within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Tests;
model Test_PEMElectrolyzerL2_CompressedStorage "Test of PEM Electrolyzer L2 connection to compressed storage, cooling model and parameters calibration"

  extends TransiEnt.Basics.Icons.Checkmodel;
  import TransiEnt;
  import Modelica.Units.SI;

  parameter String weather_file = Modelica.Utilities.Files.loadResource("modelica://H2Microgrid_TransiEnt/Resources/weather/USA_CA_Los.Angeles.Intl.AP.722950_TMY3.mos") "Path to weather file";

  // Medium declaration
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_coolant=simCenter.fluid1;

  TransiEnt.Components.Boundaries.Electrical.ActivePower.Frequency ElectricGrid_0thOrder annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-64,0})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{44,78},{64,98}})));
  StorageSystem.H2StorageSystem_Compressed H2StorageSystem(
    start_pressure=true,
    includeHeatTransfer=false,
    eta_n=0.851,
    V_geo=1.6,
    p_out=3000000,
    p_start=2950000)
                   annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={24,-60})));
  inner TransiEnt.ModelStatistics modelStatistics annotation (Placement(transformation(extent={{76,78},{96,98}})));

  Modelica.Blocks.Sources.Ramp MassflowRamp(
    offset=0,
    startTime=20,
    duration=60,
    height=1e-5) annotation (Placement(transformation(extent={{62,-38},{82,-18}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_Txim_flow H2massSink(medium=medium, variable_m_flow=true) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={58,-60})));

  Electrolyzer.Systems.SystemElectrolyzerL2_CompressedStorage ElectrolyzerSystem(
    usePowerPort=true,
    medium=medium,
    medium_coolant=medium_coolant,
    m_flow_start=1e-4,
    p_out=3000000,
    useHeatPort=true,
    useFluidCoolantPort=false,
    T_out_coolant_target=323.15,
    externalMassFlowControl=false,
    electrolyzer(voltage(humidity_const=21)),
    electrolyzer(massFlow(eta_F=1)),
    electrolyzer(pressure(p_mem_grad=17.1e5))) annotation (Placement(transformation(extent={{-14,-14},{14,14}})));
  Modelica.Blocks.Sources.CombiTimeTable Pdata(
    tableOnFile=true,
    tableName="Pdata",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/Pdata.txt"),
    verboseRead=true,
    columns={2,3,4},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=0.2) "Includes Pdata, Udata, Idata;  stop time is 4674 (timestep is 0.2 ms)"
                                                                     annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=180,
        origin={72,28})));
  Modelica.Blocks.Sources.CombiTimeTable Pstat(
    tableOnFile=true,
    tableName="Pstat",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/Pstat.txt"),
    verboseRead=true,
    columns={2,3,4},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=0.2) "Includes Pstat, Ustat, Istat; T_op_start must be 50 degC and stop time is  2280 sec (timestep is 0.2 ms)"                         annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=180,
        origin={70,56})));
  Modelica.Blocks.Sources.CombiTimeTable TempPressure(
    tableOnFile=true,
    tableName="TPdata",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/TPdata.txt"),
    verboseRead=true,
    columns={2,3},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=1) "Includes experimental H2 output Temperature and Stack Pressure over Pdata;  stop time is 23372" annotation (Placement(transformation(extent={{-90,-90},{-70,-70}})));
  Modelica.Blocks.Sources.CombiTimeTable StairSignal(table=[0,0; 499,0; 500,500; 999,500; 1000,1000; 1499,1000; 1500,1500; 1999,1500; 2000,2000; 2499,2000; 2500,2500; 2999,2500; 3000,3000; 3499,3000; 3500,3500; 3999,3500; 4000,4000; 4499,4000; 4500,4500; 4999,4500; 5000,5000; 5500,5000], tableOnFile=false) "create a stair-step signal for efficiency computation" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={42,28})));
  Modelica.Blocks.Sources.Ramp PowerRampTest(
    offset=0,
    startTime=200,
    duration=1000,
    height=5e3) "Random power curve - use P_el_set = P_el_tot" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={42,56})));
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
    relHum=0,
    TDewPoi(displayUnit="K"),
    filNam=weather_file,
    pAtmSou=Buildings.BoundaryConditions.Types.DataSource.File,
    calTSky=Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.HorizontalRadiation,
    computeWetBulbTemperature=false)
    annotation (Placement(transformation(extent={{-88,40},{-68,60}})));
  Buildings.BoundaryConditions.WeatherData.Bus
      weaBus "Weather data bus" annotation (Placement(transformation(extent={{-66,36},{-40,64}}),
                                 iconTransformation(extent={{-112,56},{-88,82}})));
equation
  connect(H2StorageSystem.H2PortOut, H2massSink.gasPort) annotation (Line(
      points={{33.8,-60},{48,-60}},
      color={255,255,0},
      thickness=1.5));
  connect(MassflowRamp.y, H2massSink.m_flow) annotation (Line(points={{83,-28},{88,-28},{88,-66},{70,-66}}, color={0,0,127}));
  connect(ElectrolyzerSystem.epp, ElectricGrid_0thOrder.epp) annotation (Line(
      points={{-14,0},{-54,0}},
      color={0,135,135},
      thickness=0.5));
  connect(ElectrolyzerSystem.gasPortOut, H2StorageSystem.H2PortIn) annotation (Line(
      points={{0,-13.86},{0,-60.1},{13.9,-60.1}},
      color={255,255,0},
      thickness=1.5));
  connect(H2StorageSystem.P_comp, ElectrolyzerSystem.CompressorPower) annotation (Line(
      points={{20.4,-70.4},{20.4,-74},{-24,-74},{-24,-8.54},{-13.58,-8.54}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(PowerRampTest.y, ElectrolyzerSystem.P_el_set) annotation (Line(points={{31,56},{0,56},{0,14.56}}, color={0,0,127}));
  connect(weaDat.weaBus, weaBus) annotation (Line(
      points={{-68,50},{-53,50}},
      color={255,204,51},
      thickness=0.5));
  connect(weaBus.TDryBul, ElectrolyzerSystem.T_environment) annotation (Line(
      points={{-52.935,50.07},{11.62,50.07},{11.62,14.42}},
      color={255,204,51},
      thickness=0.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=4600,
      Interval=1,
      __Dymola_Algorithm="Dassl"));
end Test_PEMElectrolyzerL2_CompressedStorage;