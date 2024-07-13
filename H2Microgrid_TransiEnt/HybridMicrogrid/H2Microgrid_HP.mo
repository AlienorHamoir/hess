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
                            annotation (Placement(transformation(extent={{22,-64},{62,-24}})));
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
    annotation (Placement(transformation(extent={{-58,54},{-30,80}})));
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
  Modelica.Blocks.Math.Sum PowerTotal(nin=4) annotation (Placement(transformation(extent={{-6,-8},{10,8}})));
  Modelica.Blocks.Interfaces.RealOutput SOC_battery "State of Charge battery [-]" annotation (Placement(transformation(extent={{100,72},{120,92}}), iconTransformation(extent={{100,72},{120,92}})));
  Modelica.Blocks.Interfaces.RealOutput P_battery "Battery AC power consumption [W]" annotation (Placement(transformation(extent={{100,44},{120,64}}),iconTransformation(extent={{100,44},{120,64}})));
  Modelica.Blocks.Interfaces.RealOutput P_HESS "HESS DC power production  [W]" annotation (Placement(transformation(extent={{100,-42},{120,-22}}),iconTransformation(extent={{100,-42},{120,-22}})));
  Modelica.Blocks.Interfaces.RealOutput SOC_HESS "State of Charge H2 tank  [-]" annotation (Placement(transformation(extent={{100,-62},{120,-42}}), iconTransformation(extent={{100,-62},{120,-42}})));
  Modelica.Blocks.Interfaces.RealInput P_set_HESS(unit="W", start=0) annotation (Placement(transformation(
        origin={-10,-102},
        extent={{10,-10},{-10,10}},
        rotation=270), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-106,-40})));
  Modelica.Blocks.Interfaces.RealOutput PowerBalanceMicrogrid "Net power consumption output of the microgrid = P_consumed - P_produced"
                                                                                                annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={110,0}),  iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={110,0})));
  Modelica.Blocks.Math.Gain pv_inv(k=-1)
    annotation (Placement(transformation(extent={{-4,-4},{4,4}},
        rotation=-90,
        origin={-24,38})));
  Buildings.BoundaryConditions.WeatherData.Bus
      weaBus "Weather data bus" annotation (Placement(transformation(extent={{-44,-52},{-18,-24}}),
                                 iconTransformation(extent={{-14,-70},{10,-44}})));
  Modelica.Blocks.Interfaces.RealOutput P_load "Power consumed by the load" annotation (Placement(transformation(extent={{100,-92},{120,-72}})));
  Modelica.Blocks.Interfaces.RealOutput P_PV "Power produced by the PV" annotation (Placement(transformation(extent={{100,18},{120,38}})));
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
        origin={-74,-56})));
  Modelica.Blocks.Math.Gain kWtoW(k=1000) annotation (Placement(transformation(
        extent={{-4,-4},{4,4}},
        rotation=0,
        origin={-40,-56})));
equation
  connect(pv.weaBus,weaDat. weaBus) annotation (Line(
      points={{-58,72.2},{-58,72},{-70,72}},
      color={255,204,51},
      thickness=0.5));
  connect(ctrl_PV.y,pv. scale) annotation (Line(points={{-79.6,52},{-68,52},{-68,62},{-60,62},{-60,61.8},{-60.8,61.8}},
                    color={0,0,127}));
  connect(P_set_battery,battery. PCtrl)
    annotation (Line(points={{0,106},{0,54},{16,54}},
                                                 color={0,0,127}));
  connect(battery.SOC, SOC_battery) annotation (Line(points={{62,70},{62,82},{110,82}}, color={0,0,127}));
  connect(battery.P, P_battery) annotation (Line(points={{62,54},{110,54}},                 color={0,0,127}));
  connect(PowerTotal.y, PowerBalanceMicrogrid) annotation (Line(points={{10.8,0},{110,0}},                color={0,0,127}));
  connect(hess.socTankH2, SOC_HESS) annotation (Line(points={{62.8,-50.4},{88,-50.4},{88,-52},{110,-52}}, color={0,0,127}));
  connect(pv.P, pv_inv.u) annotation (Line(points={{-28.6,67},{-24,67},{-24,42.8}}, color={0,0,127}));
  connect(P_set_HESS,hess. P_set_HESS) annotation (Line(points={{-10,-102},{-10,-50},{21,-50}},   color={0,0,127}));
  connect(hess.P_HESS, P_HESS) annotation (Line(
      points={{63.2,-31.6},{63.2,-32},{110,-32}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(pv_inv.y, PowerTotal.u[1]) annotation (Line(points={{-24,33.6},{-24,-0.6},{-7.6,-0.6}},   color={0,0,127}));
  connect(battery.P, PowerTotal.u[2]) annotation (Line(points={{62,54},{70,54},{70,22},{-14,22},{-14,0},{-10,0},{-10,-0.2},{-7.6,-0.2}},
                                                                                                                           color={0,0,127}));
  connect(weaBus, weaDat.weaBus) annotation (Line(
      points={{-31,-38},{-66,-38},{-66,72},{-70,72}},
      color={255,204,51},
      thickness=0.5));
  connect(weaBus.TDryBul,hess. T_environment) annotation (Line(
      points={{-30.935,-37.93},{-4,-37.93},{-4,-38},{21.2,-38}},
      color={255,204,51},
      thickness=0.5));
  connect(hess.P_HESS, PowerTotal.u[4]) annotation (Line(
      points={{63.2,-31.6},{72,-31.6},{72,-12},{-12,-12},{-12,0.6},{-7.6,0.6}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(pv_inv.y, P_PV) annotation (Line(points={{-24,33.6},{-24,28},{110,28}},                 color={0,0,127}));
  connect(Load.y[1], kWtoW.u) annotation (Line(points={{-63,-56},{-44.8,-56}}, color={0,0,127}));
  connect(kWtoW.y, PowerTotal.u[3]) annotation (Line(points={{-35.6,-56},{-18,-56},{-18,0.2},{-7.6,0.2}}, color={0,0,127}));
  connect(kWtoW.y, P_load) annotation (Line(points={{-35.6,-56},{-18,-56},{-18,-82},{110,-82}}, color={0,0,127}));
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
