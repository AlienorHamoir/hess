within H2Microgrid_TransiEnt;
model H2Microgrid_nonCompressed

  parameter Real building_scale = 1 "Building scale";
  parameter Real der_scale = 1 "DER scale";
  parameter Real battery_scale = 1 "Bettery scale";
  parameter Real building_ft2 = 50e3 "Building ft2 scale";
  parameter String weather_file = "" "Path to weather file";
  // parameter String weather_file = Modelica.Utilities.Files.loadResource("modelica://Buildings/Resources/weatherdata/USA_CA_San.Francisco.Intl.AP.724940_TMY3.mos") "Path to weather file";

  SCooDER.Components.Battery.Model.Battery
                                   battery(
    EMax=battery_scale*5*building_ft2*der_scale,
    Pmax=battery_scale*5/2*building_ft2*der_scale,
    SOC_start=0.1)
    annotation (Placement(transformation(extent={{20,34},{60,74}})));
  Modelica.Blocks.Interfaces.RealInput P_set_battery(unit="W", start=0)
    annotation (Placement(transformation(
        origin={0,106},
        extent={{10,-10},{-10,10}},
        rotation=90),  iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-104,52})));
  Modelica.Blocks.Math.Sum PowerTotal(nin=5) annotation (Placement(transformation(extent={{-8,6},{8,22}})));
  Modelica.Blocks.Interfaces.RealOutput SOCbattery "State of Charge battery [-]" annotation (Placement(transformation(extent={{94,72},{114,92}}), iconTransformation(extent={{94,72},{114,92}})));
  Modelica.Blocks.Interfaces.RealOutput P_battery "Battery AC power consumption [W]" annotation (Placement(transformation(extent={{94,44},{114,64}}), iconTransformation(extent={{94,44},{114,64}})));
  Modelica.Blocks.Interfaces.RealOutput P_electrolyzer "Electrolyzer and storage compressor AC power consumption [W]" annotation (Placement(transformation(extent={{94,-86},{114,-66}}), iconTransformation(extent={{94,-86},{114,-66}})));
  Modelica.Blocks.Interfaces.RealOutput P_fuelCell "Fuel cell AC power production  [W]" annotation (Placement(transformation(extent={{94,-34},{114,-14}}),iconTransformation(extent={{94,-34},{114,-14}})));
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
        origin={-106,22}), iconTransformation(
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
        origin={110,14}), iconTransformation(
        extent={{-16,-16},{16,16}},
        rotation=0,
        origin={110,14})));
  Modelica.Blocks.Math.Gain pv_inv(k=-1)
    annotation (Placement(transformation(extent={{-4,-4},{4,4}},
        rotation=-90,
        origin={-24,36})));
  HESS_nonCompressed hESS_nonCompressed annotation (Placement(transformation(extent={{20,-60},{60,-20}})));
  Modelica.Blocks.Sources.CombiTimeTable GI_PV(
    tableOnFile=true,
    tableName="PV_GI",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/PV_GI.txt"),
    verboseRead=true,
    columns={2},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=1) "GI data for PV EPFL" annotation (Placement(transformation(extent={{-100,68},{-80,88}})));
  Modelica.Blocks.Sources.CombiTimeTable Load(
    tableOnFile=true,
    tableName="Load",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/Load.txt"),
    verboseRead=true,
    columns={2},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=1) "load profiles EPFL - adapt timestep" annotation (Placement(transformation(extent={{-74,-44},{-54,-24}})));
  Modelica.Blocks.Math.Add3 PV_sum annotation (Placement(transformation(extent={{-44,62},{-32,74}})));
  Modelica.Blocks.Math.Gain Ppeak13kW(k=13000/1000) "P_peak = 13 kW for roof measurement - we use GI" annotation (Placement(transformation(extent={{-68,80},{-56,92}})));
  Modelica.Blocks.Math.Gain Ppeak16kW(k=16000/1000) "P_peak = 16 kW for roof measurement - we use GI" annotation (Placement(transformation(extent={{-68,62},{-56,74}})));
  Modelica.Blocks.Math.Gain Ppeak12kW(k=12000/1000) "P_peak = 12 kW for facade measurement - we use GVI" annotation (Placement(transformation(extent={{-68,42},{-56,54}})));
  Modelica.Blocks.Sources.CombiTimeTable GVI_PV(
    tableOnFile=true,
    tableName="PV_GVI",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/PV_GVI.txt"),
    verboseRead=true,
    columns={2},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=1) "GVI data for vertical PV EPFL" annotation (Placement(transformation(extent={{-100,38},{-80,58}})));
equation
  connect(P_set_battery,battery. PCtrl)
    annotation (Line(points={{0,106},{0,54},{16,54}},
                                                 color={0,0,127}));
  connect(battery.SOC, SOCbattery) annotation (Line(points={{62,70},{62,82},{104,82}},         color={0,0,127}));
  connect(battery.P, P_battery) annotation (Line(points={{62,54},{104,54}},                 color={0,0,127}));
  connect(PowerTotal.y, PowerBalanceMicrogrid) annotation (Line(points={{8.8,14},{110,14}},               color={0,0,127}));
  connect(pv_inv.y, P_PV) annotation (Line(points={{-24,31.6},{-24,22},{-106,22}}, color={0,0,127}));
  connect(P_set_fuelCell, hESS_nonCompressed.P_set_FC) annotation (Line(points={{-10,-102},{-10,-27.2},{19.2,-27.2}}, color={0,0,127}));
  connect(P_set_electrolyzer, hESS_nonCompressed.P_set_electrolyzer) annotation (Line(points={{8,-102},{8,-52},{18,-52}}, color={0,0,127}));
  connect(hESS_nonCompressed.P_FC, P_fuelCell) annotation (Line(
      points={{61.2,-23.6},{62,-24},{104,-24}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(hESS_nonCompressed.socTankH2, SOCtankH2) annotation (Line(points={{61.2,-40},{88,-40},{88,-48},{104,-48}}, color={0,0,127}));
  connect(hESS_nonCompressed.P_electrolyzer, P_electrolyzer) annotation (Line(
      points={{61.6,-53.6},{88,-53.6},{88,-76},{104,-76}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(hESS_nonCompressed.P_electrolyzer, PowerTotal.u[4]) annotation (Line(
      points={{61.6,-53.6},{70,-53.6},{70,-74},{-20,-74},{-20,14},{-14,14},{-14,14.32},{-9.6,14.32}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(hESS_nonCompressed.P_FC, PowerTotal.u[5]) annotation (Line(
      points={{61.2,-23.6},{68,-23.6},{68,-16},{-14,-16},{-14,14.64},{-9.6,14.64}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(Load.y[1], P_load) annotation (Line(points={{-53,-34},{-50,-34},{-50,-64},{-104,-64}}, color={0,0,127}));
  connect(Load.y[1], PowerTotal.u[3]) annotation (Line(points={{-53,-34},{-30,-34},{-30,14},{-9.6,14}}, color={0,0,127}));
  connect(pv_inv.y, PowerTotal.u[1]) annotation (Line(points={{-24,31.6},{-24,13.36},{-9.6,13.36}}, color={0,0,127}));
  connect(battery.P, PowerTotal.u[2]) annotation (Line(points={{62,54},{74,54},{74,30},{-16,30},{-16,14},{-10,14},{-10,13.68},{-9.6,13.68}}, color={0,0,127}));
  connect(GI_PV.y[1], Ppeak13kW.u) annotation (Line(points={{-79,78},{-74,78},{-74,86},{-69.2,86}}, color={0,0,127}));
  connect(GI_PV.y[1], Ppeak16kW.u) annotation (Line(points={{-79,78},{-74,78},{-74,68},{-69.2,68}}, color={0,0,127}));
  connect(GVI_PV.y[1], Ppeak12kW.u) annotation (Line(points={{-79,48},{-69.2,48}}, color={0,0,127}));
  connect(Ppeak13kW.y, PV_sum.u1) annotation (Line(points={{-55.4,86},{-52,86},{-52,72.8},{-45.2,72.8}}, color={0,0,127}));
  connect(Ppeak16kW.y, PV_sum.u2) annotation (Line(points={{-55.4,68},{-45.2,68}}, color={0,0,127}));
  connect(Ppeak12kW.y, PV_sum.u3) annotation (Line(points={{-55.4,48},{-52,48},{-52,63.2},{-45.2,63.2}}, color={0,0,127}));
  connect(PV_sum.y, pv_inv.u) annotation (Line(points={{-31.4,68},{-24,68},{-24,40.8}}, color={0,0,127}));
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
end H2Microgrid_nonCompressed;
