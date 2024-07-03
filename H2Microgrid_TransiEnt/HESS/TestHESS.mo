within H2Microgrid_TransiEnt.HESS;
model TestHESS "Testing of HESS models for microgrid applications."

  extends TransiEnt.Basics.Icons.Checkmodel;
  parameter String weather_file = Modelica.Utilities.Files.loadResource("modelica://H2Microgrid_TransiEnt/Resources/weather/USA_CA_Los.Angeles.Intl.AP.722950_TMY3.mos") "Path to weather file";


  Modelica.Blocks.Sources.Ramp PowerRampTest1(
    offset=500,
    startTime=1000,
    duration=1000,
    height=3000)
                "Random power curve - use P_el_set = P_el_tot" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-78,12})));
  HESS.HESS_Compressed hESS_Compressed annotation (Placement(transformation(extent={{6,-36},{70,18}})));
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
    relHum=0,
    TDewPoi(displayUnit="K"),
    filNam=weather_file,
    pAtmSou=Buildings.BoundaryConditions.Types.DataSource.File,
    calTSky=Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.HorizontalRadiation,
    computeWetBulbTemperature=false)
    annotation (Placement(transformation(extent={{-54,62},{-34,82}})));
  Buildings.BoundaryConditions.WeatherData.Bus
      weaBus "Weather data bus" annotation (Placement(transformation(extent={{-36,30},{-10,58}}),
                                 iconTransformation(extent={{-112,56},{-88,82}})));
  Modelica.Blocks.Sources.Sine load1(
    amplitude=3000,
    f(displayUnit="Hz") = 1,
    offset=1000)
    annotation (Placement(transformation(extent={{-86,-48},{-66,-28}})));
equation
  connect(weaBus, weaDat.weaBus) annotation (Line(
      points={{-23,44},{-22,44},{-22,72},{-34,72}},
      color={255,204,51},
      thickness=0.5));
  connect(weaBus.TDryBul, hESS_Compressed.T_environment) annotation (Line(
      points={{-22.935,44.07},{-22,44.07},{-22,-2},{2.8,-2},{2.8,-0.9}},
      color={255,204,51},
      thickness=0.5));
  connect(load1.y, hESS_Compressed.P_set_HESS) annotation (Line(points={{-65,-38},{-8,-38},{-8,-17.1},{3.12,-17.1}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>Testing of compressed and non-compressed HESS models for microgrid applications.</p>
</html>"),
    experiment(
      StopTime=1000,
      Interval=1,
      __Dymola_Algorithm="Dassl"));
end TestHESS;
