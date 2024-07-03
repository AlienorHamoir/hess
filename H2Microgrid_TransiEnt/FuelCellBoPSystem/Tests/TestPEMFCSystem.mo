within H2Microgrid_TransiEnt.FuelCellBoPSystem.Tests;
model TestPEMFCSystem "Example of a fuel cell system with its cooling and air compressor systems, with total power consumption following power setpoints"


  extends TransiEnt.Basics.Icons.Checkmodel;
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{-90,80},{-70,100}})));

  parameter String weather_file = Modelica.Utilities.Files.loadResource("modelica://H2Microgrid_TransiEnt/Resources/weather/USA_CA_Los.Angeles.Intl.AP.722950_TMY3.mos") "Path to weather file";
  parameter TransiEnt.Basics.Media.Gases.Gas_MoistAir Air=TransiEnt.Basics.Media.Gases.Gas_MoistAir() "Medium model of air" annotation (choicesAllMatching);

//     xi_const={0.01,0.7},
        //     xi_const={1,0},

  FuelCell.SystemPEMFC systemPEMFC(lambdaOController_PID(m_flow_rampup=1e-8))
                                   annotation (Placement(transformation(extent={{-26,-32},{26,20}})));
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
    relHum=0,
    TDewPoi(displayUnit="K"),
    filNam=weather_file,
    pAtmSou=Buildings.BoundaryConditions.Types.DataSource.File,
    calTSky=Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.HorizontalRadiation,
    computeWetBulbTemperature=false)
    annotation (Placement(transformation(extent={{-20,48},{0,68}})));
  Buildings.BoundaryConditions.WeatherData.Bus
      weaBus "Weather data bus" annotation (Placement(transformation(extent={{8,44},{34,72}}),
                                 iconTransformation(extent={{-112,56},{-88,82}})));
  Modelica.Blocks.Sources.Ramp PowerRamp(
    height=4000,
    duration=1000,
    offset=500,
    startTime=10) annotation (Placement(transformation(extent={{-88,-44},{-68,-24}})));
  Modelica.Blocks.Sources.Step PowerStep(
    height=1000,
    offset=0,
    startTime=1000) annotation (Placement(transformation(extent={{-88,-6},{-68,14}})));
  Modelica.Thermal.FluidHeatFlow.Examples.Utilities.DoubleRamp Load(
    startTime=1.25e3,
    interval=1.25e3,
    duration_1=1,
    offset=10,
    height_1=4000,
    height_2=-3500,
    duration_2=1) annotation (Placement(transformation(extent={{-90,-82},{-70,-62}})));
  Modelica.Blocks.Sources.Constant PowerSet(k=600) annotation (Placement(transformation(extent={{-88,34},{-68,54}})));
equation

  connect(weaDat.weaBus, weaBus) annotation (Line(
      points={{0,58},{21,58}},
      color={255,204,51},
      thickness=0.5));
  connect(systemPEMFC.T_environment, weaBus.TDryBul) annotation (Line(points={{22.88,22.08},{22.88,40},{21.065,40},{21.065,58.07}},
                                                                                                                        color={0,0,127}));
  connect(Load.y, systemPEMFC.P_el_set) annotation (Line(points={{-69,-72},{-30,-72},{-30,32},{0,32},{0,22.08}}, color={0,0,127}));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
    experiment(
      StopTime=2000,
      Interval=1,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"),
    __Dymola_experimentSetupOutput,
    __Dymola_experimentFlags(
      Advanced(
        EvaluateAlsoTop=false,
        GenerateVariableDependencies=false,
        OutputModelicaCode=false),
      Evaluate=true,
      OutputCPUtime=true,
      OutputFlatModelica=false),
    Documentation(info="<html>
<h4><span style=\"color: #008000\">1. Purpose of model</span></h4>
<p>Test environment for the PEMFC system model with different load power setpoints, accouting for temperature environment</p>
<h4><span style=\"color: #008000\">2. Level of detail, physical effects considered, and physical insight</span></h4>
<p>(Purely technical component without physical modeling.)</p>
<h4><span style=\"color: #008000\">3. Limits of validity </span></h4>
<p>(Purely technical component without physical modeling.)</p>
<h4><span style=\"color: #008000\">4. Interfaces</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">5. Nomenclature</span></h4>
<p>(no elements)</p>
<h4><span style=\"color: #008000\">6. Governing Equations</span></h4>
<p>(no equations)</p>
<h4><span style=\"color: #008000\">7. Remarks for Usage</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">8. Validation</span></h4>
<p>(no validation or testing necessary)</p>
<h4><span style=\"color: #008000\">9. References</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">10. Version History</span></h4>
</html>"));
end TestPEMFCSystem;
