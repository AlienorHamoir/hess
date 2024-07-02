within H2Microgrid_TransiEnt;
model TestHESS

    extends TransiEnt.Basics.Icons.Checkmodel;

  Modelica.Blocks.Sources.Ramp PowerRampTest1(
    offset=500,
    startTime=1000,
    duration=1000,
    height=3000)
                "Random power curve - use P_el_set = P_el_tot" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-78,12})));
  HESS_Compressed hESS_Compressed annotation (Placement(transformation(extent={{6,-36},{70,18}})));
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
    relHum=0,
    TDewPoi(displayUnit="K"),
    filNam=weather_file,
    pAtmSou=Buildings.BoundaryConditions.Types.DataSource.File,
    calTSky=Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.HorizontalRadiation,
    computeWetBulbTemperature=false)
    annotation (Placement(transformation(extent={{-30,62},{-10,82}})));
  Buildings.BoundaryConditions.WeatherData.Bus
      weaBus "Weather data bus" annotation (Placement(transformation(extent={{-12,30},{14,58}}),
                                 iconTransformation(extent={{-112,56},{-88,82}})));
  Modelica.Blocks.Sources.Sine load1(
    amplitude=3000,
    f(displayUnit="Hz") = 1,
    offset=1000)
    annotation (Placement(transformation(extent={{-86,-48},{-66,-28}})));
equation
  connect(weaBus, weaDat.weaBus) annotation (Line(
      points={{1,44},{2,44},{2,72},{-10,72}},
      color={255,204,51},
      thickness=0.5));
  connect(weaBus.TDryBul, hESS_Compressed.T_environment) annotation (Line(
      points={{1.065,44.07},{2.8,44.07},{2.8,-0.9}},
      color={255,204,51},
      thickness=0.5));
  connect(load1.y, hESS_Compressed.P_set_HESS) annotation (Line(points={{-65,-38},{-8,-38},{-8,-17.1},{3.12,-17.1}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end TestHESS;
