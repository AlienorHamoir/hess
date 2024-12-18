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
        origin={-76,-6})));
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
    amplitude=2000,
    f(displayUnit="Hz") = 0.01,
    offset=1000)
    annotation (Placement(transformation(extent={{-86,-48},{-66,-28}})));
  Modelica.Blocks.Sources.Constant StateSet(k=2)   annotation (Placement(transformation(extent={{-86,-84},{-66,-64}})));
  Modelica.Blocks.Sources.Constant PowerSet(k=0)    annotation (Placement(transformation(extent={{-88,22},{-68,42}})));
  Modelica.Blocks.Sources.Constant PowerSet1(k=0.89*4700)
                                                    annotation (Placement(transformation(extent={{-88,54},{-68,74}})));
  Modelica.Blocks.Sources.Constant StateSet1(k=0)  annotation (Placement(transformation(extent={{-32,-82},{-12,-62}})));
equation
  connect(weaBus, weaDat.weaBus) annotation (Line(
      points={{-23,44},{-22,44},{-22,72},{-34,72}},
      color={255,204,51},
      thickness=0.5));
  connect(weaBus.TDryBul, hESS_Compressed.T_environment) annotation (Line(
      points={{-22.935,44.07},{-22,44.07},{-22,-2},{6,-2},{6,-9}},
      color={255,204,51},
      thickness=0.5));
  connect(StateSet.y,hESS_Compressed.state_FC)  annotation (Line(points={{-65,-74},{-40,-74},{-40,0},{-8,0},{-8,4.5},{6,4.5}},
                                                                                                          color={0,0,127}));
  connect(PowerSet.y, hESS_Compressed.P_set_EL) annotation (Line(points={{-67,32},{-42,32},{-42,-24.66},{6.32,-24.66}},
                                                                                                                      color={0,0,127}));
  connect(PowerSet1.y, hESS_Compressed.P_set_FC) annotation (Line(points={{-67,64},{-58,64},{-58,86},{-6,86},{-6,13.14},{6.32,13.14}},
                                                                                                                                   color={0,0,127}));
  connect(StateSet1.y, hESS_Compressed.state_EL) annotation (Line(points={{-11,-72},{6,-72},{6,-32.76}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>Testing of compressed and non-compressed HESS models for microgrid applications.</p>
</html>"),
    experiment(
      StopTime=1000,
      Interval=1,
      __Dymola_Algorithm="Dassl"));
end TestHESS;
