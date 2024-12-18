within H2Microgrid_TransiEnt.HybridMicrogrid;
model H2Microgrid_HP "Hybrid microgrid with high-pressure compressed storage"

  parameter Real building_scale = 1 "Building scale";
  parameter Real der_scale = 0.25 "DER scale - assumption in our case, we have 1 floor and PV on 25% of roof surface";
  parameter Real battery_scale = 1 "Battery scale";
  parameter Real building_ft2 = 5500 "Building ft2 scale";
  parameter Real sqft2sqm = 10.765 "Convert ft^2 surface to m^2";

  parameter Modelica.Units.SI.Pressure p_max=350e5   "Maximum pressure in storage" annotation (Dialog(tab="General", group="H2 storage parameters"));
  parameter Real LOH_start_HESS=0.5  "SOC of H2 tank storage at t=0" annotation (Dialog(tab="General", group="Initial parameters"));
  parameter Real SOC_start_battery=0.5  "SOC of battery at t=0" annotation (Dialog(tab="General", group="Initial parameters"));


  parameter String weather_file = Modelica.Utilities.Files.loadResource("modelica://H2Microgrid_TransiEnt/Resources/weather/USA_CA_Los.Angeles.Intl.AP.722950_TMY3.mos") "Path to weather file";
  parameter String load_file = ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/loads/commercial_SmallOffice_LA.txt") "Path to load file";
  parameter String disturbance_file = ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/loads/LoadDisturbance_Year.txt") "Path to load disturbance file";

  HESS.HESS_Compressed hess(p_max=p_max, SOC_start=LOH_start_HESS)
                            annotation (Placement(transformation(extent={{20,-44},{60,-4}})));
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
    relHum=0,
    TDewPoi(displayUnit="K"),
    filNam=weather_file,
    pAtmSou=Buildings.BoundaryConditions.Types.DataSource.File,
    calTSky=Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.HorizontalRadiation,
    computeWetBulbTemperature=false)
    annotation (Placement(transformation(extent={{-90,44},{-70,64}})));
  SCooDER.Components.Photovoltaics.Model.PVModule_simple
                                                 pv(n=der_scale*building_ft2/(1.65*sqft2sqm))
    annotation (Placement(transformation(extent={{-58,36},{-34,60}})));
  Modelica.Blocks.Sources.Constant ctrl_PV(k=1)
    annotation (Placement(transformation(extent={{-88,30},{-80,38}})));
  SCooDER.Components.Battery.Model.Battery
                                   battery(
    EMax=5000,
    Pmax=5000,
    SOC_start=SOC_start_battery,
    SOC_min=0.1,
    SOC_max=0.9,
    etaCha=0.98,
    etaDis=0.98)
    annotation (Placement(transformation(extent={{20,30},{60,70}})));
  Modelica.Blocks.Interfaces.RealInput P_set_battery(unit="W", start=0)
    annotation (Placement(transformation(
        origin={0,80},
        extent={{10,-10},{-10,10}},
        rotation=90),  iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-104,52})));
  Modelica.Blocks.Math.Sum GridPowerTotal(nin=5) annotation (Placement(transformation(extent={{2,2},{18,18}})));
  Modelica.Blocks.Interfaces.RealOutput P_battery "Battery AC power consumption [W]" annotation (Placement(transformation(extent={{100,40},{120,60}}),iconTransformation(extent={{100,40},{120,60}})));
  Modelica.Blocks.Interfaces.RealOutput P_FC "Fuel cell system DC power balance  [W]" annotation (Placement(transformation(extent={{100,-14},{120,6}}),   iconTransformation(extent={{100,-14},{120,6}})));
  Modelica.Blocks.Interfaces.RealOutput LOH "H2 tank Level of Hydrogen [-]" annotation (Placement(transformation(extent={{100,-30},{120,-10}}), iconTransformation(extent={{100,-30},{120,-10}})));
  Modelica.Blocks.Interfaces.RealInput P_set_FC(unit="W", start=0) annotation (Placement(transformation(
        origin={-16,-86},
        extent={{10,-10},{-10,10}},
        rotation=270), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-106,-30})));
  Modelica.Blocks.Math.Gain pv_inv(k=-1)
    annotation (Placement(transformation(extent={{-4,-4},{4,4}},
        rotation=-90,
        origin={-24,20})));
  Buildings.BoundaryConditions.WeatherData.Bus
      weaBus "Weather data bus" annotation (Placement(transformation(extent={{-52,-38},{-26,-10}}),
                                 iconTransformation(extent={{-14,-70},{10,-44}})));
  Modelica.Blocks.Interfaces.RealOutput P_load "Power consumed by the load" annotation (Placement(transformation(extent={{100,-68},{120,-48}}), iconTransformation(extent={{100,-68},{120,-48}})));
  Modelica.Blocks.Interfaces.RealOutput P_PV "Power produced by the PV" annotation (Placement(transformation(extent={{100,18},{120,38}}), iconTransformation(extent={{100,18},{120,38}})));
  Modelica.Blocks.Math.Gain kWtoW(k=1000) annotation (Placement(transformation(
        extent={{-4,-4},{4,4}},
        rotation=0,
        origin={-30,-58})));
  Modelica.Blocks.Interfaces.RealInput P_set_EL(unit="W", start=0) annotation (Placement(transformation(
        origin={6,-86},
        extent={{10,-10},{-10,10}},
        rotation=270), iconTransformation(
        extent={{10,10},{-10,-10}},
        rotation=180,
        origin={-106,-70})));
  Modelica.Blocks.Interfaces.RealOutput P_EL "Electrolyzer system DC power balance  [W]" annotation (Placement(transformation(extent={{100,-48},{120,-28}}), iconTransformation(extent={{100,-48},{120,-28}})));
  Modelica.Blocks.Interfaces.RealInput state_FC "FC state computed by controller" annotation (Placement(transformation(
        extent={{-12,-12},{12,12}},
        rotation=90,
        origin={-50,-86}),  iconTransformation(
        extent={{-12,-12},{12,12}},
        rotation=90,
        origin={-50,-86})));
  Modelica.Blocks.Interfaces.RealInput state_EL "Electrolyzer state computed by controller" annotation (Placement(transformation(
        extent={{-12,-12},{12,12}},
        rotation=90,
        origin={48,-88}),  iconTransformation(
        extent={{-12,-12},{12,12}},
        rotation=90,
        origin={48,-88})));
  Modelica.Blocks.Interfaces.RealOutput SOE_battery "State of Energy [Wh]" annotation (Placement(transformation(extent={{100,58},{120,78}})));
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
        origin={-80,-40})));
  Modelica.Blocks.Sources.CombiTimeTable Disturbance(
    tableOnFile=true,
    tableName="Load",
    fileName=disturbance_file,
    verboseRead=true,
    columns={2},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=60) "Base on load file from Matthieu " annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-80,-70})));
  Modelica.Blocks.Math.Add add annotation (Placement(transformation(extent={{-62,-64},{-48,-50}})));
equation
  connect(pv.weaBus,weaDat. weaBus) annotation (Line(
      points={{-58,52.8},{-58,54},{-70,54}},
      color={255,204,51},
      thickness=0.5));
  connect(ctrl_PV.y,pv. scale) annotation (Line(points={{-79.6,34},{-68,34},{-68,44},{-60,44},{-60,43.2},{-60.4,43.2}},
                    color={0,0,127}));
  connect(P_set_battery,battery. PCtrl)
    annotation (Line(points={{0,80},{0,50},{16,50}},
                                                 color={0,0,127}));
  connect(battery.P, P_battery) annotation (Line(points={{62,50},{110,50}},                 color={0,0,127}));
  connect(hess.LOH, LOH) annotation (Line(points={{60.8,-21.2},{86,-21.2},{86,-20},{110,-20}},
                                                                         color={0,0,127}));
  connect(pv.P, pv_inv.u) annotation (Line(points={{-32.8,48},{-24,48},{-24,24.8}}, color={0,0,127}));
  connect(P_set_FC, hess.P_set_FC) annotation (Line(points={{-16,-86},{-16,-7.6},{20.2,-7.6}},color={0,0,127}));
  connect(hess.P_FC, P_FC) annotation (Line(
      points={{61,-11.4},{61,-4},{110,-4}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(pv_inv.y, GridPowerTotal.u[1]) annotation (Line(points={{-24,15.6},{-24,9.36},{0.4,9.36}},   color={0,0,127}));
  connect(battery.P, GridPowerTotal.u[2]) annotation (Line(points={{62,50},{86,50},{86,24},{-14,24},{-14,8},{-10,8},{-10,9.68},{0.4,9.68}},   color={0,0,127}));
  connect(weaBus, weaDat.weaBus) annotation (Line(
      points={{-39,-24},{-66,-24},{-66,54},{-70,54}},
      color={255,204,51},
      thickness=0.5));
  connect(weaBus.TDryBul,hess. T_environment) annotation (Line(
      points={{-38.935,-23.93},{-4,-23.93},{-4,-24},{20,-24}},
      color={255,204,51},
      thickness=0.5));
  connect(kWtoW.y, GridPowerTotal.u[3]) annotation (Line(points={{-25.6,-58},{-20,-58},{-20,10},{0.4,10}},
                                                                                                         color={0,0,127}));
  connect(kWtoW.y, P_load) annotation (Line(points={{-25.6,-58},{110,-58}},                     color={0,0,127}));
  connect(P_set_EL, hess.P_set_EL) annotation (Line(points={{6,-86},{6,-35.6},{20.2,-35.6}},
                                                                                          color={0,0,127}));
  connect(hess.P_EL, P_EL) annotation (Line(
      points={{61.4,-37},{94,-37},{94,-38},{110,-38}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(hess.P_FC, GridPowerTotal.u[4]) annotation (Line(
      points={{61,-11.4},{70,-11.4},{70,2},{-8,2},{-8,10.32},{0.4,10.32}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(hess.P_EL, GridPowerTotal.u[5]) annotation (Line(
      points={{61.4,-37},{70,-37},{70,-52},{-20,-52},{-20,10.64},{0.4,10.64}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(state_FC,hess.state_FC)  annotation (Line(points={{-50,-86},{-50,-70},{10,-70},{10,-54},{20,-54},{20,-14}},  color={0,0,127}));
  connect(state_EL,hess.state_EL)  annotation (Line(points={{48,-88},{48,-41.6},{20,-41.6}},
                                                                                 color={0,0,127}));
  connect(battery.SOE, SOE_battery) annotation (Line(points={{62,60},{84,60},{84,68},{110,68}}, color={0,0,127}));
  connect(pv.P, P_PV) annotation (Line(points={{-32.8,48},{-24,48},{-24,28},{110,28}},                                 color={0,0,127}));
  connect(Load.y[1], add.u1) annotation (Line(points={{-69,-40},{-68,-40},{-68,-46},{-63.4,-46},{-63.4,-52.8}}, color={0,0,127}));
  connect(Disturbance.y[1], add.u2) annotation (Line(points={{-69,-70},{-63.4,-70},{-63.4,-61.2}}, color={0,0,127}));
  connect(add.y, kWtoW.u) annotation (Line(points={{-47.3,-57},{-47.3,-58},{-34.8,-58}}, color={0,0,127}));
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
HP")}),                            Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<h4>1. Purpose of model</h4>
<p>Hybrid microgrid, containing a HESS with compressed storage, a BESS, PV generation and buildings loads inputs.</p>
<p>Meant to be later connected to a MPC as an FMU.</p>
<h4>2. Level of detail, physical effects considered, and physical insight</h4>
<p>PV and BESS models from SCooDER library. The BESS model is a simplified battery model.</p>
<p>Parameters and data are coming from the DOE or from the DESL setup.</p>
<h4>3. Limits of validity </h4>
<h4>4. Interfaces</h4>
<h4>5. Nomenclature</h4>
<p>(no remarks)</p>
<h4>6. Governing Equations</h4>
<h4>7. Remarks for Usage</h4>
<h4>8. Validation</h4>
<p>Tested in &quot;H2Microgrid_TransiEnt.HybridMicrogrid.TestMicrogrid&quot;</p>
<h4>9. References</h4>
<h4>10. Version History</h4>
<p>Model created by Ali&eacute;nor Hamoir in June 2024</p>
</html>"),
    experiment(
      StopTime=40000000,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"));
end H2Microgrid_HP;
