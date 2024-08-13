within H2Microgrid_TransiEnt.FuelCellBoPSystem.Tests;
model TestPEMFCSystemexp "Example of a fuel cell system with its cooling and air compressor systems, with total power consumption following power setpoints"

  extends TransiEnt.Basics.Icons.Checkmodel;

  parameter String weather_file = Modelica.Utilities.Files.loadResource("modelica://H2Microgrid_TransiEnt/Resources/weather/USA_CA_Los.Angeles.Intl.AP.722950_TMY3.mos") "Path to weather file";
  parameter TransiEnt.Basics.Media.Gases.Gas_MoistAir Air=TransiEnt.Basics.Media.Gases.Gas_MoistAir() "Medium model of air" annotation (choicesAllMatching);

//     xi_const={0.01,0.7},
        //     xi_const={1,0},

  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
    relHum=0,
    TDewPoi(displayUnit="K"),
    filNam=weather_file,
    pAtmSou=Buildings.BoundaryConditions.Types.DataSource.File,
    calTSky=Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.HorizontalRadiation,
    computeWetBulbTemperature=false)
    annotation (Placement(transformation(extent={{-40,72},{-20,92}})));
  Buildings.BoundaryConditions.WeatherData.Bus
      weaBus "Weather data bus" annotation (Placement(transformation(extent={{-16,68},{10,96}}),
                                 iconTransformation(extent={{-112,56},{-88,82}})));
  Modelica.Blocks.Sources.Ramp PowerRamp(
    height=4700,
    duration=15,
    offset=300,
    startTime=1)  annotation (Placement(transformation(extent={{-88,-44},{-68,-24}})));
  Modelica.Blocks.Sources.Step PowerStep(
    height=650,
    offset=0,
    startTime=1000) "can go up to 680 W with curret parameters"
                    annotation (Placement(transformation(extent={{-88,-6},{-68,14}})));
  Modelica.Thermal.FluidHeatFlow.Examples.Utilities.DoubleRamp Load(
    startTime=1.25e3,
    interval=1.25e3,
    duration_1=1000,
    offset=0,
    height_1=5000,
    height_2=-4700,
    duration_2=1000)
                  annotation (Placement(transformation(extent={{-88,-80},{-68,-60}})));
  Modelica.Blocks.Sources.Constant PowerSet(k=-600)
                                                   annotation (Placement(transformation(extent={{-88,34},{-68,54}})));
  Modelica.Blocks.Sources.CombiTimeTable StairSignal(table=[0,0; 499,0; 500,500; 999,500; 1000,1000; 1499,1000; 1500,1500; 1999,1500; 2000,2000; 2499,2000; 2500,2500; 2999,2500; 3000,3000; 3499,3000; 3500,3500; 3999,3500; 4000,4000; 4499,4000; 4500,4500; 4999,4500; 5000,5000; 5500,5000], tableOnFile=false) "create a stair-step signal for efficiency computation" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={76,78})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{-88,76},{-68,96}})));
  Modelica.Blocks.Sources.CombiTimeTable Pstat(
    tableOnFile=true,
    tableName="Pstat",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/Pstat.txt"),
    verboseRead=true,
    columns={2,3,4},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=0.2) "Includes Pstat, Ustat, Istat; T_op_start must be 50 degC and stop time is  2280 sec (timestep is 0.2 ms)"                         annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=180,
        origin={76,46})));
  FuelCell.SystemPEMFCexp systemPEMFCexp(
    coolingModel(k_p=0.5),
    powerController(i_max=300),
    FC(T_stack(start=333.15))) annotation (Placement(transformation(extent={{-40,-20},{0,20}})));
  Modelica.Blocks.Math.Gain A_cell(k=0.2) annotation (Placement(transformation(
        extent={{-6,-6},{6,6}},
        rotation=180,
        origin={48,78})));
  Modelica.Blocks.Math.Add add annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={38,40})));
  Modelica.Blocks.Sources.Constant const(k=-2900) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={76,16})));
  Modelica.Blocks.Math.Gain A_cell1(k=0.5)
                                          annotation (Placement(transformation(
        extent={{-6,-6},{6,6}},
        rotation=180,
        origin={8,40})));
  Modelica.Blocks.Sources.Sine HESScommand(
    amplitude=1000,
    f(displayUnit="Hz") = 0.0001,
    offset=1500) annotation (Placement(transformation(extent={{-54,-80},{-34,-60}})));
equation

  connect(weaDat.weaBus, weaBus) annotation (Line(
      points={{-20,82},{-3,82}},
      color={255,204,51},
      thickness=0.5));
  connect(systemPEMFCexp.T_environment, weaBus.TDryBul) annotation (Line(points={{-2.4,21.6},{-2,21.6},{-2,82.07},{-2.935,82.07}}, color={0,0,127}), Text(
      string="%second",
      index=1,
      extent={{-3,6},{-3,6}},
      horizontalAlignment=TextAlignment.Right));
  connect(StairSignal.y[1], A_cell.u) annotation (Line(points={{65,78},{55.2,78}}, color={0,0,127}));
  connect(Pstat.y[1], add.u2) annotation (Line(points={{65,46},{50,46}}, color={0,0,127}));
  connect(const.y, add.u1) annotation (Line(points={{65,16},{50,16},{50,34}}, color={0,0,127}));
  connect(add.y, A_cell1.u) annotation (Line(points={{27,40},{15.2,40}}, color={0,0,127}));
  connect(PowerRamp.y, systemPEMFCexp.P_el_set) annotation (Line(points={{-67,-34},{-44,-34},{-44,32},{-20,32},{-20,21.6}}, color={0,0,127}));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
    experiment(
      StopTime=6500,
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
end TestPEMFCSystemexp;
