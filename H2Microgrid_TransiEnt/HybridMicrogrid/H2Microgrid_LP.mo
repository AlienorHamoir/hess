within H2Microgrid_TransiEnt.HybridMicrogrid;
model H2Microgrid_LP "Hybrid microgrid with low-pressure non-compressed storage"

  parameter Real building_scale = 1 "Building scale";
  parameter Real der_scale = 1 "DER scale";
  parameter Real battery_scale = 1 "Bettery scale";
  parameter Real building_ft2 = 50e3 "Building ft2 scale";
  parameter String weather_file = "" "Path to weather file";
  // parameter String weather_file = Modelica.Utilities.Files.loadResource("modelica://Buildings/Resources/weatherdata/USA_CA_San.Francisco.Intl.AP.724940_TMY3.mos") "Path to weather file";

  Modelica.Blocks.Interfaces.RealInput P_set_battery(unit="W", start=0)
    annotation (Placement(transformation(
        origin={0,106},
        extent={{10,-10},{-10,10}},
        rotation=90),  iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-104,52})));
  Modelica.Blocks.Math.Sum PowerTotal(nin=4) annotation (Placement(transformation(extent={{-6,-8},{10,8}})));
  Modelica.Blocks.Interfaces.RealOutput SOCbattery "State of Charge battery [-]" annotation (Placement(transformation(extent={{94,72},{114,92}}), iconTransformation(extent={{94,72},{114,92}})));
  Modelica.Blocks.Interfaces.RealOutput P_battery "Battery AC power consumption [W]" annotation (Placement(transformation(extent={{96,44},{116,64}}), iconTransformation(extent={{96,44},{116,64}})));
  Modelica.Blocks.Interfaces.RealOutput P_HESS "HESS AC power balance  [W]" annotation (Placement(transformation(extent={{94,-40},{114,-20}}), iconTransformation(extent={{94,-40},{114,-20}})));
  Modelica.Blocks.Interfaces.RealOutput SOCtankH2 "State of Charge H2 tank  [-]" annotation (Placement(transformation(extent={{94,-64},{114,-44}}), iconTransformation(extent={{94,-64},{114,-44}})));
  Modelica.Blocks.Interfaces.RealInput P_set_HESS(unit="W", start=0) annotation (Placement(transformation(
        origin={-10,-102},
        extent={{10,-10},{-10,10}},
        rotation=270), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-106,-40})));
  Modelica.Blocks.Interfaces.RealOutput PowerBalanceMicrogrid "Net power consumption output of the microgrid = P_consumed - P_produced"
                                                                                                annotation (Placement(transformation(
        extent={{-16,-16},{16,16}},
        rotation=0,
        origin={110,6}),  iconTransformation(
        extent={{-16,-16},{16,16}},
        rotation=0,
        origin={110,6})));
  Modelica.Blocks.Math.Gain pv_inv(k=-1)
    annotation (Placement(transformation(extent={{-4,-4},{4,4}},
        rotation=-90,
        origin={-24,38})));
  HESS.HESS_nonCompressed hESS_nonCompressed annotation (Placement(transformation(extent={{20,-66},{60,-26}})));
  Modelica.Blocks.Sources.CombiTimeTable GI_PV(
    tableOnFile=true,
    tableName="PV_GI",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/PV_GI.txt"),
    verboseRead=true,
    columns={2},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=1) "GI data for PV EPFL" annotation (Placement(transformation(extent={{-100,62},{-80,82}})));
  Modelica.Blocks.Math.Add3 PV_sum annotation (Placement(transformation(extent={{-44,56},{-32,68}})));
  Modelica.Blocks.Math.Gain Ppeak13kW(k=13000/1000) "P_peak = 13 kW for roof measurement - we use GI" annotation (Placement(transformation(extent={{-68,74},{-56,86}})));
  Modelica.Blocks.Math.Gain Ppeak16kW(k=16000/1000) "P_peak = 16 kW for roof measurement - we use GI" annotation (Placement(transformation(extent={{-68,56},{-56,68}})));
  Modelica.Blocks.Math.Gain Ppeak12kW(k=12000/1000) "P_peak = 12 kW for facade measurement - we use GVI" annotation (Placement(transformation(extent={{-68,36},{-56,48}})));
  Modelica.Blocks.Sources.CombiTimeTable GVI_PV(
    tableOnFile=true,
    tableName="PV_GVI",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/PV_GVI.txt"),
    verboseRead=true,
    columns={2},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=1) "GVI data for vertical PV EPFL" annotation (Placement(transformation(extent={{-100,32},{-80,52}})));
  Modelica.Blocks.Sources.CombiTimeTable Load(
    tableOnFile=true,
    tableName="Load",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/Load.txt"),
    verboseRead=true,
    columns={2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=1) "Includes 21 different building load profiles; unit is kW? ;  stop time is 28800" annotation (Placement(transformation(extent={{-84,-50},{-64,-30}})));
  Modelica.Blocks.Math.Gain kWtoW(k=1000) "Load profiles are given in kW, we exploit the microgrid in W" annotation (Placement(transformation(extent={{-52,-46},{-40,-34}})));
  SCooDER.Components.Battery.Model.Battery
                                   battery(
    EMax=25000,
    Pmax=25000,
    SOC_start=0.5,
    SOC_min=0.1,
    SOC_max=0.9,
    etaCha=0.98,
    etaDis=0.98)
    annotation (Placement(transformation(extent={{28,46},{68,86}})));
  Modelica.Blocks.Sources.CombiTimeTable Load_DOE(
    tableOnFile=true,
    tableName="Load",
    fileName=load_file,
    verboseRead=true,
    columns={2},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=1) "Load profiles from DOE" annotation (Placement(transformation(extent={{-84,-18},{-64,2}})));
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
    relHum=0,
    TDewPoi(displayUnit="K"),
    filNam=weather_file,
    pAtmSou=Buildings.BoundaryConditions.Types.DataSource.File,
    calTSky=Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.HorizontalRadiation,
    computeWetBulbTemperature=false)
    annotation (Placement(transformation(extent={{-94,-82},{-74,-62}})));
  Buildings.BoundaryConditions.WeatherData.Bus
      weaBus "Weather data bus" annotation (Placement(transformation(extent={{-56,-86},{-30,-58}}),
                                 iconTransformation(extent={{-14,-70},{10,-44}})));
equation
  connect(PowerTotal.y, PowerBalanceMicrogrid) annotation (Line(points={{10.8,0},{60,0},{60,6},{110,6}},  color={0,0,127}));
  connect(P_set_HESS, hESS_nonCompressed.P_set_HESS) annotation (Line(points={{-10,-102},{-10,-52.2},{17.4,-52.2}}, color={0,0,127}));
  connect(hESS_nonCompressed.P_HESS, P_HESS) annotation (Line(
      points={{61.6,-30.8},{61.6,-30},{104,-30}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(hESS_nonCompressed.socTankH2, SOCtankH2) annotation (Line(points={{61.8,-46.2},{88,-46.2},{88,-54},{104,-54}},
                                                                                                                     color={0,0,127}));
  connect(GI_PV.y[1], Ppeak13kW.u) annotation (Line(points={{-79,72},{-74,72},{-74,80},{-69.2,80}}, color={0,0,127}));
  connect(GI_PV.y[1], Ppeak16kW.u) annotation (Line(points={{-79,72},{-74,72},{-74,62},{-69.2,62}}, color={0,0,127}));
  connect(GVI_PV.y[1], Ppeak12kW.u) annotation (Line(points={{-79,42},{-69.2,42}}, color={0,0,127}));
  connect(Ppeak13kW.y, PV_sum.u1) annotation (Line(points={{-55.4,80},{-52,80},{-52,66.8},{-45.2,66.8}}, color={0,0,127}));
  connect(Ppeak16kW.y, PV_sum.u2) annotation (Line(points={{-55.4,62},{-45.2,62}}, color={0,0,127}));
  connect(Ppeak12kW.y, PV_sum.u3) annotation (Line(points={{-55.4,42},{-52,42},{-52,57.2},{-45.2,57.2}}, color={0,0,127}));
  connect(PV_sum.y, pv_inv.u) annotation (Line(points={{-31.4,62},{-24,62},{-24,42.8}}, color={0,0,127}));
  connect(Load.y[1], kWtoW.u) annotation (Line(points={{-63,-40},{-53.2,-40}}, color={0,0,127}));
  connect(P_set_battery, battery.PCtrl) annotation (Line(points={{0,106},{0,66},{24,66}}, color={0,0,127}));
  connect(battery.SOC, SOCbattery) annotation (Line(points={{70,82},{104,82}}, color={0,0,127}));
  connect(battery.P, P_battery) annotation (Line(points={{70,66},{92,66},{92,54},{106,54}}, color={0,0,127}));
  connect(pv_inv.y, PowerTotal.u[1]) annotation (Line(points={{-24,33.6},{-24,-0.6},{-7.6,-0.6}}, color={0,0,127}));
  connect(battery.P, PowerTotal.u[2]) annotation (Line(points={{70,66},{78,66},{78,28},{-16,28},{-16,0},{-12,0},{-12,-0.2},{-7.6,-0.2}}, color={0,0,127}));
  connect(kWtoW.y, PowerTotal.u[3]) annotation (Line(points={{-39.4,-40},{-20,-40},{-20,0.2},{-7.6,0.2}}, color={0,0,127}));
  connect(hESS_nonCompressed.P_HESS, PowerTotal.u[4]) annotation (Line(
      points={{61.6,-30.8},{62,-30.8},{62,-14},{-12,-14},{-12,0.6},{-7.6,0.6}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(weaDat.weaBus, weaBus) annotation (Line(
      points={{-74,-72},{-43,-72}},
      color={255,204,51},
      thickness=0.5), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}},
      horizontalAlignment=TextAlignment.Left));
  connect(weaBus.TDryBul, hESS_nonCompressed.T_environment) annotation (Line(
      points={{-42.935,-71.93},{-12,-71.93},{-12,-41.6},{17.2,-41.6}},
      color={255,204,51},
      thickness=0.5), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}},
      horizontalAlignment=TextAlignment.Right));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,128,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid), Text(
          extent={{-78,44},{76,-42}},
          textColor={102,44,145},
          fontName="Arial Black",
          textStyle={TextStyle.Bold},
          textString="H2 Grid
LP")}),                            Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<h4>1. Purpose of model</h4>
<p>Hybrid microgrid, containing a HESS with non-compressed storage, a BESS, PV generation and buildings loads inputs.</p>
<p>Meant to be later connected to a MPC as an FMU.</p>
<h4>2. Level of detail, physical effects considered, and physical insight</h4>
<p>BESS model is a simplified battery model from SCooDER library. </p>
<p>Parameters and data are coming from the DOE or from the DESL setup.</p>
<h4>3. Limits of validity </h4>
<h4>4. Interfaces</h4>
<h4>5. Nomenclature</h4>
<h4>6. Governing Equations</h4>
<h4>7. Remarks for Usage</h4>
<h4>8. Validation</h4>
<p>Tested in &quot;H2Microgrid_TransiEnt.HybridMicrogrid.TestMicrogrid&quot;</p>
<h4>9. References</h4>
<h4>10. Version History</h4>
<p>Model created by Ali&eacute;nor Hamoir in June 2024</p>
</html>"));
end H2Microgrid_LP;
