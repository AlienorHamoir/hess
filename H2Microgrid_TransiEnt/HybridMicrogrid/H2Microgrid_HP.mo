within H2Microgrid_TransiEnt.HybridMicrogrid;
model H2Microgrid_HP "Hybrid microgrid with high-pressure compressed storage"

  parameter Real building_scale = 1 "Building scale";
  parameter Real der_scale = 0.25 "DER scale - assumption in our case, we have 1 floor and PV on 25% of roof surface";
  parameter Real battery_scale = 1 "Battery scale";
  parameter Real building_ft2 = 5500 "Building ft2 scale";
  parameter Real sqft2sqm = 10.765 "Convert ft^2 surface to m^2";

  parameter Modelica.Units.SI.Pressure p_max=350e5   "Maximum pressure in storage" annotation (Dialog(tab="General", group="H2 storage parameters"));
  parameter Real SOC_start_HESS=0.5  "SOC of H2 tank storage at t=0" annotation (Dialog(tab="General", group="Initial parameters"));
  parameter Real SOC_start_battery=0.5  "SOC of battery at t=0" annotation (Dialog(tab="General", group="Initial parameters"));


  parameter String weather_file = Modelica.Utilities.Files.loadResource("modelica://H2Microgrid_TransiEnt/Resources/weather/USA_CA_Los.Angeles.Intl.AP.722950_TMY3.mos") "Path to weather file";
  parameter String load_file = ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/loads/commercial_SmallOffice_LA.txt") "Path to load file";

  HESS.HESS_Compressed hess(p_max=p_max, SOC_start=SOC_start_HESS)
                            annotation (Placement(transformation(extent={{20,-60},{60,-20}})));
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
    relHum=0,
    TDewPoi(displayUnit="K"),
    filNam=weather_file,
    pAtmSou=Buildings.BoundaryConditions.Types.DataSource.File,
    calTSky=Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.HorizontalRadiation,
    computeWetBulbTemperature=false)
    annotation (Placement(transformation(extent={{-90,62},{-70,82}})));
  SCooDER.Components.Photovoltaics.Model.PVModule_simple
                                                 pv(n=der_scale*building_ft2/(1.65*sqft2sqm))
    annotation (Placement(transformation(extent={{-58,54},{-34,78}})));
  Modelica.Blocks.Sources.Constant ctrl_PV(k=1)
    annotation (Placement(transformation(extent={{-88,48},{-80,56}})));
  SCooDER.Components.Battery.Model.Battery
                                   battery(
    EMax=5000,
    Pmax=5000,
    SOC_start=SOC_start_battery,
    SOC_min=0.1,
    SOC_max=0.9,
    etaCha=0.98,
    etaDis=0.98)
    annotation (Placement(transformation(extent={{20,34},{60,74}})));
  Modelica.Blocks.Interfaces.RealInput P_set_battery(unit="W", start=0)
    annotation (Placement(transformation(
        origin={0,106},
        extent={{10,-10},{-10,10}},
        rotation=90),  iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-104,52})));
  Modelica.Blocks.Math.Sum PowerTotal(nin=5) annotation (Placement(transformation(extent={{0,-8},{16,8}})));
  Modelica.Blocks.Interfaces.RealOutput SOC_battery "State of Charge battery [-]" annotation (Placement(transformation(extent={{100,80},{120,100}}),iconTransformation(extent={{100,80},{120,100}})));
  Modelica.Blocks.Interfaces.RealOutput P_battery "Battery AC power consumption [W]" annotation (Placement(transformation(extent={{100,60},{120,80}}),iconTransformation(extent={{100,60},{120,80}})));
  Modelica.Blocks.Interfaces.RealOutput P_FC "Fuel cell system DC power balance  [W]" annotation (Placement(transformation(extent={{100,-30},{120,-10}}), iconTransformation(extent={{100,-30},{120,-10}})));
  Modelica.Blocks.Interfaces.RealOutput SOC_HESS "State of Charge H2 tank  [-]" annotation (Placement(transformation(extent={{100,-50},{120,-30}}), iconTransformation(extent={{100,-50},{120,-30}})));
  Modelica.Blocks.Interfaces.RealInput P_set_FC(unit="W", start=0) annotation (Placement(transformation(
        origin={-16,-102},
        extent={{10,-10},{-10,10}},
        rotation=270), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-106,-30})));
  Modelica.Blocks.Interfaces.RealOutput PowerBalanceMicrogrid "Net power consumption output of the microgrid = P_consumed - P_produced"
                                                                                                annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={110,10}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={110,10})));
  Modelica.Blocks.Math.Gain pv_inv(k=-1)
    annotation (Placement(transformation(extent={{-4,-4},{4,4}},
        rotation=-90,
        origin={-24,38})));
  Buildings.BoundaryConditions.WeatherData.Bus
      weaBus "Weather data bus" annotation (Placement(transformation(extent={{-52,-54},{-26,-26}}),
                                 iconTransformation(extent={{-14,-70},{10,-44}})));
  Modelica.Blocks.Interfaces.RealOutput P_load "Power consumed by the load" annotation (Placement(transformation(extent={{100,-98},{120,-78}}), iconTransformation(extent={{100,-98},{120,-78}})));
  Modelica.Blocks.Interfaces.RealOutput P_PV "Power produced by the PV" annotation (Placement(transformation(extent={{100,40},{120,60}}), iconTransformation(extent={{100,40},{120,60}})));
  Modelica.Blocks.Sources.CombiTimeTable Load(
    tableOnFile=true,
    tableName="Load",
    fileName=load_file,
    verboseRead=true,
    columns={2},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=60) "Base load in LA, from DOE"                                                                                annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-84,-74})));
  Modelica.Blocks.Math.Gain kWtoW(k=1000) annotation (Placement(transformation(
        extent={{-4,-4},{4,4}},
        rotation=0,
        origin={-50,-74})));
  Modelica.Blocks.Interfaces.RealInput P_set_EL(unit="W", start=0) annotation (Placement(transformation(
        origin={6,-102},
        extent={{10,-10},{-10,10}},
        rotation=270), iconTransformation(
        extent={{10,10},{-10,-10}},
        rotation=180,
        origin={-106,-70})));
  Modelica.Blocks.Interfaces.RealOutput P_EL "Electrolyzer system DC power balance  [W]" annotation (Placement(transformation(extent={{100,-74},{120,-54}}), iconTransformation(extent={{100,-70},{120,-50}})));
  Modelica.Blocks.Interfaces.RealInput state_FC "FC state computed by controller" annotation (Placement(transformation(
        extent={{-12,-12},{12,12}},
        rotation=90,
        origin={-50,-102}), iconTransformation(
        extent={{-12,-12},{12,12}},
        rotation=90,
        origin={-50,-102})));
  Modelica.Blocks.Interfaces.RealInput state_EL "Electrolyzer state computed by controller" annotation (Placement(transformation(
        extent={{-12,-12},{12,12}},
        rotation=90,
        origin={48,-104}), iconTransformation(
        extent={{-12,-12},{12,12}},
        rotation=90,
        origin={48,-104})));
equation
  connect(pv.weaBus,weaDat. weaBus) annotation (Line(
      points={{-58,70.8},{-58,72},{-70,72}},
      color={255,204,51},
      thickness=0.5));
  connect(ctrl_PV.y,pv. scale) annotation (Line(points={{-79.6,52},{-68,52},{-68,62},{-60,62},{-60,61.2},{-60.4,61.2}},
                    color={0,0,127}));
  connect(P_set_battery,battery. PCtrl)
    annotation (Line(points={{0,106},{0,54},{16,54}},
                                                 color={0,0,127}));
  connect(battery.SOC, SOC_battery) annotation (Line(points={{62,70},{62,90},{110,90}}, color={0,0,127}));
  connect(battery.P, P_battery) annotation (Line(points={{62,54},{86,54},{86,70},{110,70}}, color={0,0,127}));
  connect(PowerTotal.y, PowerBalanceMicrogrid) annotation (Line(points={{16.8,0},{64,0},{64,10},{110,10}},color={0,0,127}));
  connect(hess.SOC, SOC_HESS) annotation (Line(points={{60.8,-40},{110,-40}}, color={0,0,127}));
  connect(pv.P, pv_inv.u) annotation (Line(points={{-32.8,66},{-24,66},{-24,42.8}}, color={0,0,127}));
  connect(P_set_FC, hess.P_set_FC) annotation (Line(points={{-16,-102},{-16,-26},{19.8,-26}}, color={0,0,127}));
  connect(hess.P_FC, P_FC) annotation (Line(
      points={{61.8,-26.6},{61.8,-20},{110,-20}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(pv_inv.y, PowerTotal.u[1]) annotation (Line(points={{-24,33.6},{-24,-0.64},{-1.6,-0.64}}, color={0,0,127}));
  connect(battery.P, PowerTotal.u[2]) annotation (Line(points={{62,54},{70,54},{70,22},{-14,22},{-14,0},{-10,0},{-10,-0.32},{-1.6,-0.32}},
                                                                                                                           color={0,0,127}));
  connect(weaBus, weaDat.weaBus) annotation (Line(
      points={{-39,-40},{-66,-40},{-66,72},{-70,72}},
      color={255,204,51},
      thickness=0.5));
  connect(weaBus.TDryBul,hess. T_environment) annotation (Line(
      points={{-38.935,-39.93},{-4,-39.93},{-4,-40},{20,-40}},
      color={255,204,51},
      thickness=0.5));
  connect(pv_inv.y, P_PV) annotation (Line(points={{-24,33.6},{-24,50},{110,50}},                 color={0,0,127}));
  connect(Load.y[1], kWtoW.u) annotation (Line(points={{-73,-74},{-54.8,-74}}, color={0,0,127}));
  connect(kWtoW.y, PowerTotal.u[3]) annotation (Line(points={{-45.6,-74},{-20,-74},{-20,0},{-1.6,0}},     color={0,0,127}));
  connect(kWtoW.y, P_load) annotation (Line(points={{-45.6,-74},{94,-74},{94,-88},{110,-88}},   color={0,0,127}));
  connect(P_set_EL, hess.P_set_EL) annotation (Line(points={{6,-102},{6,-54},{20.2,-54}}, color={0,0,127}));
  connect(hess.P_EL, P_EL) annotation (Line(
      points={{61.8,-53.8},{94,-53.8},{94,-64},{110,-64}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(hess.P_FC, PowerTotal.u[4]) annotation (Line(
      points={{61.8,-26.6},{70,-26.6},{70,-14},{-8,-14},{-8,0.32},{-1.6,0.32}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(hess.P_EL, PowerTotal.u[5]) annotation (Line(
      points={{61.8,-53.8},{70,-53.8},{70,-68},{-20,-68},{-20,0.64},{-1.6,0.64}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(state_FC, hess.delta_FC) annotation (Line(points={{-50,-102},{-50,-86},{10,-86},{10,-70},{32,-70},{32,-60}}, color={0,0,127}));
  connect(state_EL, hess.delta_EL) annotation (Line(points={{48,-104},{48,-60}}, color={0,0,127}));
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
    experiment(StopTime=40000000, __Dymola_Algorithm="Dassl"));
end H2Microgrid_HP;
