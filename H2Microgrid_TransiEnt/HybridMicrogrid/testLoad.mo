within H2Microgrid_TransiEnt.HybridMicrogrid;
model testLoad

      extends TransiEnt.Basics.Icons.Checkmodel;

  parameter Real building_scale = 1 "Building scale";
  parameter Real der_scale = 0.25 "DER scale - assumption in our case, we have 1 floor and PV on 50% of roof surface";
  parameter Real battery_scale = 1 "Battery scale";
  parameter Real building_ft2 = 5500 "Building ft2 scale";
  parameter Real sqft2sqm = 10.765 "Convert ft^2 surface to m^2";
  parameter String weather_file = Modelica.Utilities.Files.loadResource("modelica://H2Microgrid_TransiEnt/Resources/weather/USA_CA_Los.Angeles.Intl.AP.722950_TMY3.mos") "Path to weather file";
  parameter String load_file = ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/loads/commercial_SmallOffice_LA.txt") "Path to load file";



  Modelica.Blocks.Sources.CombiTimeTable Load(
    tableOnFile=true,
    tableName="Load",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/loads/commercial_SmallOffice_LA.txt"),
    verboseRead=true,
    columns={2},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=60) "Base load in LA, from DOE"                                                                                annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-10,-28})));
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
    relHum=0,
    TDewPoi(displayUnit="K"),
    filNam=weather_file,
    pAtmSou=Buildings.BoundaryConditions.Types.DataSource.File,
    calTSky=Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.HorizontalRadiation,
    computeWetBulbTemperature=false)
    annotation (Placement(transformation(extent={{-38,36},{-18,56}})));
  SCooDER.Components.Photovoltaics.Model.PVModule_simple
                                                 pv(n=der_scale*building_ft2/(1.65*sqft2sqm))
    annotation (Placement(transformation(extent={{-6,28},{22,54}})));
  Modelica.Blocks.Sources.Constant ctrl_PV(k=1)
    annotation (Placement(transformation(extent={{-36,22},{-28,30}})));
  Modelica.Blocks.Math.Gain kWtoW(k=1000) annotation (Placement(transformation(
        extent={{-4,-4},{4,4}},
        rotation=0,
        origin={20,-28})));
equation
  connect(pv.weaBus,weaDat. weaBus) annotation (Line(
      points={{-6,46.2},{-6,46},{-18,46}},
      color={255,204,51},
      thickness=0.5));
  connect(ctrl_PV.y, pv.scale) annotation (Line(points={{-27.6,26},{-8.8,26},{-8.8,35.8}}, color={0,0,127}));
  connect(Load.y[1], kWtoW.u) annotation (Line(points={{1,-28},{15.2,-28}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=35552000,
      __Dymola_NumberOfIntervals=1000,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"));
end testLoad;
