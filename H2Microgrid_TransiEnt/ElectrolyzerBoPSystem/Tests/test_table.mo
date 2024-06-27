within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Tests;
model test_table
  Modelica.Blocks.Sources.CombiTimeTable combiTimeTable_Pdata(
    tableOnFile=true,
    tableName="Pdata",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/Pdata.txt"),
    verboseRead=true,
    columns={2,3,4},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=1)                                                "stop time is 23000"
                                                                annotation (Placement(transformation(extent={{-20,4},{0,24}})));
  Modelica.Blocks.Sources.CombiTimeTable combiTimeTable_TempPres(
    tableOnFile=true,
    tableName="TPdata",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/TPdata.txt"),
    verboseRead=true,
    columns={2,3},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=1) "Includes Temperature and Stack pressure over Pdata;  stop time is 23372" annotation (Placement(transformation(extent={{-20,46},{0,66}})));
  Modelica.Blocks.Sources.CombiTimeTable combiTimeTable_Pstat(
    tableOnFile=true,
    tableName="Pstat",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/Pstat.txt"),
    verboseRead=true,
    columns={2,3,4},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments) "Includes Pstat, Ustat, Istat; T_start must be around 40 degC and stop time is 13680 sec"
                                                                                                 annotation (Placement(transformation(extent={{-22,-54},{-2,-34}})));
  Modelica.Blocks.Sources.CombiTimeTable combiTimeTable_StairSignal(table=[ 0, 0;
   499, 0;
   500, 500;
   999, 500;
   1000, 1000;
   1499, 1000;
   1500, 1500;
   1999, 1500;
   2000, 2000;
   2499, 2000;
   2500, 2500;
   2999, 2500;
   3000, 3000;
   3499, 3000;
   3500, 3500;
   3999, 3500;
   4000, 4000;
   4499, 4000;
   4500, 4500;
   4999, 4500;
   5000, 5000;
   5500, 5000], tableOnFile=false) "create a stair-step signal" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={32,-10})));
  Modelica.Blocks.Sources.CombiTimeTable Load(
    tableOnFile=true,
    tableName="Load",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://H2Microgrid_TransiEnt/Resources/Load.txt"),
    verboseRead=true,
    columns={2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    timeScale=1) "Includes 21 different buildig load profiles; unit is kW? ;  stop time is 28800" annotation (Placement(transformation(extent={{-78,34},{-58,54}})));
  Modelica.Blocks.Sources.Ramp ramp(
    height=-600,
    duration=100,
    offset=-200,
    startTime=10) annotation (Placement(transformation(extent={{40,20},{60,40}})));
  Modelica.Blocks.Sources.Ramp ramp1(
    height=-300,
    duration=100,
    offset=-500,
    startTime=10) annotation (Placement(transformation(extent={{76,-22},{96,-2}})));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=6000,
      Interval=1,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"));
end test_table;
