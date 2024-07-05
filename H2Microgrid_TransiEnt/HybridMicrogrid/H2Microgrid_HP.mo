within H2Microgrid_TransiEnt.HybridMicrogrid;
model H2Microgrid_HP "Hybrid microgrid with high-pressure compressed storage"

  parameter Real building_scale = 1 "Building scale";
  parameter Real der_scale = 1 "DER scale";
  parameter Real battery_scale = 1 "Bettery scale";
  parameter Real building_ft2 = 50e3 "Building ft2 scale";
  parameter String weather_file = Modelica.Utilities.Files.loadResource("modelica://H2Microgrid_TransiEnt/Resources/weather/USA_CA_Los.Angeles.Intl.AP.722950_TMY3.mos") "Path to weather file";
//   parameter String load_file = Modelica.Utilities.Files.loadResource("modelica://H2Microgrid_TransiEnt/Resources/loads/USA_CA_Los.Angeles.Intl.AP.722950_TMY3.mos") "Path to weather file";

  HESS.HESS_Compressed hESS annotation (Placement(transformation(extent={{22,-64},{62,-24}})));
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
    relHum=0,
    TDewPoi(displayUnit="K"),
    filNam=weather_file,
    pAtmSou=Buildings.BoundaryConditions.Types.DataSource.File,
    calTSky=Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.HorizontalRadiation,
    computeWetBulbTemperature=false)
    annotation (Placement(transformation(extent={{-90,62},{-70,82}})));
  SCooDER.Components.Photovoltaics.Model.PVModule_simple
                                                 pv(n=0.004*building_ft2*der_scale)
    annotation (Placement(transformation(extent={{-58,54},{-30,80}})));
  Modelica.Blocks.Sources.Constant ctrl_PV(k=1)
    annotation (Placement(transformation(extent={{-88,48},{-80,56}})));
  SCooDER.Components.Battery.Model.Battery
                                   battery(
    EMax=25000,
    Pmax=25000,
    SOC_start=0.5,
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
  Modelica.Blocks.Interfaces.RealOutput SOCbattery "State of Charge battery [-]" annotation (Placement(transformation(extent={{100,72},{120,92}}),iconTransformation(extent={{100,72},{120,92}})));
  Modelica.Blocks.Interfaces.RealOutput P_battery "Battery AC power consumption [W]" annotation (Placement(transformation(extent={{100,44},{120,64}}),iconTransformation(extent={{100,44},{120,64}})));
  Modelica.Blocks.Interfaces.RealOutput P_HESS "HESS DC power production  [W]" annotation (Placement(transformation(extent={{100,-36},{120,-16}}),iconTransformation(extent={{100,-36},{120,-16}})));
  Modelica.Blocks.Interfaces.RealOutput SOCtankH2 "State of Charge H2 tank  [-]" annotation (Placement(transformation(extent={{100,-60},{120,-40}}),iconTransformation(extent={{100,-60},{120,-40}})));
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
        origin={110,8}),  iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={110,8})));
  Modelica.Blocks.Math.Gain pv_inv(k=-1)
    annotation (Placement(transformation(extent={{-4,-4},{4,4}},
        rotation=-90,
        origin={-24,38})));
  Buildings.BoundaryConditions.WeatherData.Bus
      weaBus "Weather data bus" annotation (Placement(transformation(extent={{-38,-52},{-12,-24}}),
                                 iconTransformation(extent={{-14,-70},{10,-44}})));
  Modelica.Blocks.Sources.Constant Load(k=0) annotation (Placement(transformation(extent={{-90,-36},{-70,-16}})));
  Modelica.Blocks.Interfaces.RealOutput P_load "Power consumed by the load" annotation (Placement(transformation(extent={{100,-92},{120,-72}})));
  Modelica.Blocks.Interfaces.RealOutput P_PV "Power produced by the PV" annotation (Placement(transformation(extent={{100,22},{120,42}})));
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
  connect(battery.SOC, SOCbattery) annotation (Line(points={{62,70},{62,82},{110,82}},         color={0,0,127}));
  connect(battery.P, P_battery) annotation (Line(points={{62,54},{110,54}},                 color={0,0,127}));
  connect(PowerTotal.y, PowerBalanceMicrogrid) annotation (Line(points={{10.8,0},{58,0},{58,8},{110,8}},  color={0,0,127}));
  connect(hESS.socTankH2, SOCtankH2) annotation (Line(points={{62.8,-50.4},{88,-50.4},{88,-50},{110,-50}},
                                                                                                       color={0,0,127}));
  connect(pv.P, pv_inv.u) annotation (Line(points={{-28.6,67},{-24,67},{-24,42.8}}, color={0,0,127}));
  connect(P_set_HESS, hESS.P_set_HESS) annotation (Line(points={{-10,-102},{-10,-50},{20.2,-50}}, color={0,0,127}));
  connect(hESS.P_HESS, P_HESS) annotation (Line(
      points={{63.2,-31.6},{63.2,-30},{86,-30},{86,-26},{110,-26}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(pv_inv.y, PowerTotal.u[1]) annotation (Line(points={{-24,33.6},{-24,-0.6},{-7.6,-0.6}},   color={0,0,127}));
  connect(battery.P, PowerTotal.u[2]) annotation (Line(points={{62,54},{70,54},{70,22},{-14,22},{-14,0},{-10,0},{-10,-0.2},{-7.6,-0.2}},
                                                                                                                           color={0,0,127}));
  connect(weaBus, weaDat.weaBus) annotation (Line(
      points={{-25,-38},{-66,-38},{-66,72},{-70,72}},
      color={255,204,51},
      thickness=0.5));
  connect(weaBus.TDryBul, hESS.T_environment) annotation (Line(
      points={{-24.935,-37.93},{-4,-37.93},{-4,-38},{20,-38}},
      color={255,204,51},
      thickness=0.5));
  connect(hESS.P_HESS, PowerTotal.u[4]) annotation (Line(
      points={{63.2,-31.6},{66,-31.6},{66,-12},{-18,-12},{-18,0.6},{-7.6,0.6}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(Load.y, PowerTotal.u[3]) annotation (Line(points={{-69,-26},{-44,-26},{-44,0.2},{-7.6,0.2}}, color={0,0,127}));
  connect(Load.y, P_load) annotation (Line(points={{-69,-26},{-44,-26},{-44,-82},{110,-82}}, color={0,0,127}));
  connect(pv_inv.y, P_PV) annotation (Line(points={{-24,33.6},{-24,28},{94,28},{94,32},{110,32}}, color={0,0,127}));
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
</html>"));
end H2Microgrid_HP;
