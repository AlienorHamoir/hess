within H2Microgrid_TransiEnt.HybridMicrogrid;
model testLoad

      extends TransiEnt.Basics.Icons.Checkmodel;


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
        origin={-10,10})));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=3153600,
      Interval=1,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"));
end testLoad;
