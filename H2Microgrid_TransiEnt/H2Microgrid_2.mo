within H2Microgrid_TransiEnt;
model H2Microgrid_2

  parameter Real building_scale = 1 "Building scale";
  parameter Real der_scale = 1 "DER scale";
  parameter Real battery_scale = 1 "Bettery scale";
  parameter Real building_ft2 = 50e3 "Building ft2 scale";
  parameter String weather_file = "" "Path to weather file";
  // parameter String weather_file = Modelica.Utilities.Files.loadResource("modelica://Buildings/Resources/weatherdata/USA_CA_San.Francisco.Intl.AP.724940_TMY3.mos") "Path to weather file";


  HESS hESS annotation (Placement(transformation(extent={{20,-60},{60,-20}})));
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
    relHum=0,
    TDewPoi(displayUnit="K"),
    filNam=weather_file,
    pAtmSou=Buildings.BoundaryConditions.Types.DataSource.File,
    calTSky=Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.HorizontalRadiation,
    computeWetBulbTemperature=false)
    annotation (Placement(transformation(extent={{-92,42},{-72,62}})));
  SCooDER.Components.Photovoltaics.Model.PVModule_simple
                                                 pv(n=0.004*building_ft2*der_scale)
    annotation (Placement(transformation(extent={{-60,34},{-32,60}})));
  Modelica.Blocks.Sources.Constant ctrl_PV(k=1)
    annotation (Placement(transformation(extent={{-90,28},{-82,36}})));
  SCooDER.Components.Battery.Model.Battery
                                   battery(
    EMax=battery_scale*5*building_ft2*der_scale,
    Pmax=battery_scale*5/2*building_ft2*der_scale,
    SOC_start=0.1)
    annotation (Placement(transformation(extent={{20,20},{60,60}})));
  Modelica.Blocks.Interfaces.RealInput P_set_battery(unit="W", start=0)
    annotation (Placement(transformation(
        origin={0,106},
        extent={{10,-10},{-10,10}},
        rotation=90),  iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-104,52})));
  Modelica.Blocks.Math.Sum PowerTotal(nin=5) annotation (Placement(transformation(extent={{-8,-8},{8,8}})));
  Modelica.Blocks.Sources.Constant Load(k=0) annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));
  Modelica.Blocks.Interfaces.RealOutput SOCbattery "State of Charge battery [-]" annotation (Placement(transformation(extent={{94,72},{114,92}}), iconTransformation(extent={{94,72},{114,92}})));
  Modelica.Blocks.Interfaces.RealOutput P_battery "Battery AC power consumption [W]" annotation (Placement(transformation(extent={{94,44},{114,64}}), iconTransformation(extent={{94,44},{114,64}})));
  Modelica.Blocks.Interfaces.RealOutput P_electrolyzer "Electrolyzer and storage compressor AC power consumption [W]" annotation (Placement(transformation(extent={{94,-86},{114,-66}}), iconTransformation(extent={{94,-86},{114,-66}})));
  Modelica.Blocks.Interfaces.RealOutput P_fuelCell "Fuel cell AC power production  [W]" annotation (Placement(transformation(extent={{94,-28},{114,-8}}), iconTransformation(extent={{94,-28},{114,-8}})));
  Modelica.Blocks.Interfaces.RealOutput SOCtankH2 "State of Charge H2 tank  [-]" annotation (Placement(transformation(extent={{94,-58},{114,-38}}), iconTransformation(extent={{94,-58},{114,-38}})));
  Modelica.Blocks.Interfaces.RealOutput P_load "Load AC power consumption [W]" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-104,-64}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={70,100})));
  Modelica.Blocks.Interfaces.RealOutput P_PV "PV AC power production [W]" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-106,16}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={28,100})));
  Modelica.Blocks.Interfaces.RealInput P_set_electrolyzer(unit="W", start=0) annotation (Placement(transformation(
        origin={8,-102},
        extent={{10,-10},{-10,10}},
        rotation=270), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-106,0})));
  Modelica.Blocks.Interfaces.RealInput P_set_fuelCell(unit="W", start=0) annotation (Placement(transformation(
        origin={-14,-102},
        extent={{10,-10},{-10,10}},
        rotation=270), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-106,-40})));
  Modelica.Blocks.Interfaces.RealOutput PowerBalanceMicrogrid "Connector of Real output signal" annotation (Placement(transformation(
        extent={{-16,-16},{16,16}},
        rotation=0,
        origin={110,14}), iconTransformation(
        extent={{-16,-16},{16,16}},
        rotation=0,
        origin={110,14})));
  Modelica.Blocks.Math.Gain pv_inv(k=-1)
    annotation (Placement(transformation(extent={{-4,-4},{4,4}},
        rotation=-90,
        origin={-24,28})));
  Modelica.Blocks.Math.Gain pv_inv1(k=-1)
    annotation (Placement(transformation(extent={{-4,-4},{4,4}},
        rotation=0,
        origin={72,-18})));
equation
  connect(pv.weaBus,weaDat. weaBus) annotation (Line(
      points={{-60,52.2},{-60,52},{-72,52}},
      color={255,204,51},
      thickness=0.5));
  connect(ctrl_PV.y,pv. scale) annotation (Line(points={{-81.6,32},{-62.8,32},{-62.8,41.8}},
                    color={0,0,127}));
  connect(P_set_battery,battery. PCtrl)
    annotation (Line(points={{0,106},{0,40},{16,40}},
                                                 color={0,0,127}));
  connect(battery.P, PowerTotal.u[2]) annotation (Line(points={{62,40},{72,40},{72,16},{-18,16},{-18,-0.32},{-9.6,-0.32}}, color={0,0,127}));
  connect(Load.y, PowerTotal.u[3]) annotation (Line(points={{-59,-30},{-28,-30},{-28,0},{-9.6,0}}, color={0,0,127}));
  connect(battery.SOC, SOCbattery) annotation (Line(points={{62,56},{84,56},{84,82},{104,82}}, color={0,0,127}));
  connect(battery.P, P_battery) annotation (Line(points={{62,40},{84,40},{84,54},{104,54}}, color={0,0,127}));
  connect(Load.y, P_load) annotation (Line(points={{-59,-30},{-56,-30},{-56,-64},{-104,-64}}, color={0,0,127}));
  connect(PowerTotal.y, PowerBalanceMicrogrid) annotation (Line(points={{8.8,0},{88,0},{88,14},{110,14}}, color={0,0,127}));
  connect(P_fuelCell, PowerTotal.u[4]) annotation (Line(points={{104,-18},{-16,-18},{-16,0.32},{-9.6,0.32}},                   color={0,0,127}));
  connect(P_electrolyzer, PowerTotal.u[5]) annotation (Line(points={{104,-76},{-22,-76},{-22,0},{-9.6,0},{-9.6,0.64}}, color={0,0,127}));
  connect(P_set_electrolyzer, hESS.P_set_electrolyzer) annotation (Line(points={{8,-102},{8,-46.4},{19.6,-46.4}}, color={0,0,127}));
  connect(hESS.socTankH2, SOCtankH2) annotation (Line(points={{61.2,-40},{88,-40},{88,-48},{104,-48}}, color={0,0,127}));
  connect(pv.P, pv_inv.u) annotation (Line(points={{-30.6,47},{-24,47},{-24,32.8}}, color={0,0,127}));
  connect(pv_inv.y, PowerTotal.u[1]) annotation (Line(points={{-24,23.6},{-24,0},{-12,0},{-12,-0.64},{-9.6,-0.64}}, color={0,0,127}));
  connect(pv_inv.y, P_PV) annotation (Line(points={{-24,23.6},{-24,16},{-106,16}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,128,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid), Text(
          extent={{-78,44},{76,-42}},
          textColor={102,44,145},
          fontName="Arial Black",
          textStyle={TextStyle.Bold},
          textString="H2 Grid")}), Diagram(coordinateSystem(preserveAspectRatio=false)));
end H2Microgrid_2;
