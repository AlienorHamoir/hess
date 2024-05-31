within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Tests;
model TestPEMElectrolyzerL1_Dryer_Storage "Test of PEM Electrolyzer L1 connection to drying system and storage"

  extends TransiEnt.Basics.Icons.Checkmodel;
  import TransiEnt;
  import      Modelica.Units.SI;

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_gas=simCenter.gasModel4;

//protected
  function plotResult

  constant String resultFileName = "TestPEMElectrolyzer_L1_Dynamics.mat";

  output String resultFile;

  algorithm
    clearlog();
    assert(cd(Modelica.Utilities.System.getEnvironmentVariable(TransiEnt.Basics.Types.WORKINGDIR)), "Error changing directory: Working directory must be set as environment variable with name 'workingdir' for this script to work.");
    resultFile :=TransiEnt.Basics.Functions.fullPathName(Modelica.Utilities.System.getEnvironmentVariable(TransiEnt.Basics.Types.WORKINGDIR) + "/" + resultFileName);
    removePlots(false);
    createPlot(id=1, position={0, 0, 1563, 749}, y={"electrolyzer_0thOrder.dynamics.H_flow_H2", "electrolyzer_1stOrder.dynamics.H_flow_H2","electrolyzer_2ndOrder.dynamics.H_flow_H2"},range={0.0,52.0,-100000.0,1100000.0},erase=false,grid=true,filename=resultFile,colors={{28,108,200},{238,46,47},{0,140,72}});
  end plotResult;

public
  parameter SI.Power P_el_n=1e6 "Nominal electrical power of the electrolyzer";
  parameter SI.Power P_el_min=0.05*P_el_n "Minimal electrical power of the electrolyzer";
  parameter SI.Power P_el_max=1.0*P_el_n "Maximal electrical power of the electrolyzer";
  parameter SI.Temperature T_out=273.15+10 "Temperature of the produced hydrogen";
  parameter SI.Efficiency eta_n=0.75 "Nominal efficiency of the electrolyzer";
  parameter SI.Pressure p_out=30e5 "Pressure of the produced hydrogen";

  TransiEnt.Components.Boundaries.Electrical.ActivePower.Frequency ElectricGrid_0thOrder annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-72,22})));
  Modelica.Blocks.Sources.Ramp PowerRamp(
    offset=2,
    startTime=0,
    duration=20,
    height=40000) annotation (Placement(transformation(extent={{-94,70},{-74,90}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y=5) annotation (Placement(transformation(extent={{20,40},{40,60}})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{-88,-96},{-68,-76}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_pTxi boundary_pTxi(medium=medium_gas) annotation (Placement(transformation(extent={{82,-72},{62,-52}})));
  inner TransiEnt.ModelStatistics modelStatistics annotation (Placement(transformation(extent={{-54,-96},{-34,-76}})));
  H2DryingSystem.H2DryingSubsystem h2DryingSubsystem annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={-12,-42})));
  StorageSystem.H2StorageSystem_nonCompressed h2StorageSystem_nonCompressed annotation (Placement(transformation(extent={{12,-72},{32,-52}})));
equation
  connect(ElectricGrid_0thOrder.epp, electrolyzerSystem.epp) annotation (Line(
      points={{-62,22},{-60,21},{-34,21}},
      color={0,135,135},
      thickness=0.5));
  connect(PowerRamp.y, electrolyzerSystem.P_el_set) annotation (Line(points={{-73,80},{-11,80},{-11,42.84}},
                                                                                                           color={0,0,127}));
  connect(realExpression.y, electrolyzerSystem.m_flow_feedIn) annotation (Line(points={{41,50},{46,50},{46,64},{12,64},{12,37.8}},color={0,0,127}));
  connect(electrolyzerSystem.gasPortOut, h2DryingSubsystem.gasPortIn) annotation (Line(
      points={{-11,0.21},{-12,0.21},{-12,-32.2}},
      color={255,255,0},
      thickness=1.5));
  connect(h2DryingSubsystem.gasPortOut, h2StorageSystem_nonCompressed.H2PortIn) annotation (Line(
      points={{-12,-52},{-12,-62},{12.6,-62}},
      color={255,255,0},
      thickness=1.5));
  connect(h2StorageSystem_nonCompressed.H2PortOut, boundary_pTxi.gasPort) annotation (Line(
      points={{33,-62},{62,-62}},
      color={255,255,0},
      thickness=1.5));
  connect(electrolyzerSystem.pressureIn, h2StorageSystem_nonCompressed.pressureTank) annotation (Line(points={{2.8,44.1},{2.8,68},{48,68},{48,-56.6},{32.6,-56.6}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(StopTime=1000, __Dymola_Algorithm="Dassl"));
end TestPEMElectrolyzerL1_Dryer_Storage;
