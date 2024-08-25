within H2Microgrid_TransiEnt.FuelCellBoPSystem.Tests;
model TestPEMFCexp "Test and validation with experimental results of the PEMFC model"

//________________________________________________________________________________//
// Component of the TransiEnt Library, version: 2.0.3                             //
//                                                                                //
// Licensed by Hamburg University of Technology under the 3-BSD-clause.           //
// Copyright 2021, Hamburg University of Technology.                              //
//________________________________________________________________________________//
//                                                                                //
// TransiEnt.EE, ResiliEntEE, IntegraNet and IntegraNet II are research projects  //
// supported by the German Federal Ministry of Economics and Energy               //
// (FKZ 03ET4003, 03ET4048, 0324027 and 03EI1008).                                //
// The TransiEnt Library research team consists of the following project partners://
// Institute of Engineering Thermodynamics (Hamburg University of Technology),    //
// Institute of Energy Systems (Hamburg University of Technology),                //
// Institute of Electrical Power and Energy Technology                            //
// (Hamburg University of Technology)                                             //
// Fraunhofer Institute for Environmental, Safety, and Energy Technology UMSICHT, //
// Gas- und WÃ¤rme-Institut Essen                                                  //
// and                                                                            //
// XRG Simulation GmbH (Hamburg, Germany).                                        //
//________________________________________________________________________________//

  extends TransiEnt.Basics.Icons.Checkmodel;

 // parameter TransiEnt.Basics.Media.Gases.Gas_MoistAir Air=TransiEnt.Basics.Media.Gases.Gas_MoistAir() "Medium model of air" annotation (choicesAllMatching);

  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow AirSource(
    variable_m_flow=true,
    variable_xi=false,
    m_flow_const=0.001,
    T_const=23.5 + 273,
    medium=FC.Air) annotation (Placement(transformation(
        extent={{6.5,-9},{-6.5,9}},
        rotation=180,
        origin={-39.5,-27})));
//     xi_const={0.01,0.7},
  ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi AirSink(
    variable_p=false,
    variable_T=false,
    variable_xi=false,
    p_const=1e5,
    T_const=23.50 + 273.15,
    medium=FC.Air) annotation (Placement(transformation(
        extent={{-7,-8},{7,8}},
        rotation=180,
        origin={61,-28})));
        //     xi_const={1,0},
  ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi SyngasSink(
    variable_p=false,
    variable_T=false,
    variable_xi=false,
    p_const=1e5,
    medium=FC.Syngas,
    T_const=23.5 + 273.15,
    xi_const={0,0,0,0,1,0}) annotation (Placement(transformation(
        extent={{-6,-9},{6,9}},
        rotation=180,
        origin={60,1})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow SyngasSource(
    variable_T=false,
    m_flow_const=5.1e-2,
    variable_m_flow=true,
    variable_xi=false,
    T_const=40 + 273.15,
    medium=FC.Syngas,
    xi_const={0,0,0,0,1,0}) annotation (Placement(transformation(extent={{-48,-9},{-32,7}})));

  Modelica.Blocks.Sources.Ramp PolarizationCurrentRamp(
    height=117,
    duration=2.34,
    offset=8,
    startTime=0) "Laurencelle Ballard MK5-E polarization curve" annotation (Placement(transformation(extent={{-66,22},{-50,38}})));
  Modelica.Thermal.FluidHeatFlow.Examples.Utilities.DoubleRamp DynamicCurrentDensity(
    startTime=13,
    interval=26,
    duration_1=1,
    offset=0.015,
    height_1=0.095,
    height_2=-0.095,
    duration_2=1) "Test current density for dynamical model validation (Khan et Iqbal)"
                                                                        annotation (Placement(transformation(extent={{-34,74},{-18,90}})));
  Modelica.Blocks.Sources.Ramp CurrentRangeRamp(
    height=0.8,
    duration=200,
    offset=0.005,
    startTime=1) "Current density operating range" annotation (Placement(transformation(extent={{-92,74},{-76,90}})));
  Modelica.Blocks.Sources.Constant TempSet(k=15 + 273.15) annotation (Placement(transformation(extent={{0,28},{8,36}})));
  Modelica.Blocks.Sources.Step CurrentDensityStep(
    height=0.75,
    offset=0.015,
    startTime=100)
                  annotation (Placement(transformation(extent={{-62,74},{-46,90}})));
  FuelCell.PEMFC FC(
    usePowerPort=false,
    useHeatPort=false,
    T_stack(displayUnit="degC", start=345.15)) annotation (Placement(transformation(extent={{-6,-16},{14,4}})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{72,-94},{88,-78}})));
  Modelica.Blocks.Math.Gain A_cell(k=232) annotation (Placement(transformation(
        extent={{-5,-5},{5,5}},
        rotation=-90,
        origin={-29,49})));
  Modelica.Blocks.Sources.CombiTimeTable CellVoltageExp(table=[0,0.995; 13,0.995; 14,0.82; 16,0.8; 28,0.8; 40,0.8; 41,0.97; 42,0.993; 43,0.995; 70,0.995], tableOnFile=false) "cell dynamic voltage response (Kan et Iqbal)" annotation (Placement(transformation(
        extent={{-7,-7},{7,7}},
        rotation=180,
        origin={83,89})));
  Controller.LambdaController_PID lambdaController_PID(lambda_target=1.5, m_flow_rampup=1e-8) annotation (Placement(transformation(extent={{-90,-10},{-70,10}})));
  Controller.LambdaController_PID lambdaController_PID1(lambda_target=2,    m_flow_rampup=1e-8) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-26,-54})));
  Modelica.Thermal.FluidHeatFlow.Examples.Utilities.DoubleRamp LoadChange(
    startTime=4,
    interval=0.1,
    duration_1=0.1,
    offset=100,
    height_1=120,
    height_2=-65,
    duration_2=0.3) "Test current for dynamical model validation (Laurencelle)"
                                                                          annotation (Placement(transformation(extent={{-92,46},{-76,62}})));
  CoolingSystem.HeatPortCooling.CoolingModel coolingModel(
    k_p=0.5,                                                        tau_i=0.01,
    controller(controllerType=Modelica.Blocks.Types.SimpleController.P))        annotation (Placement(transformation(extent={{20,22},{40,42}})));
  Modelica.Blocks.Sources.CombiTimeTable StackVoltageExp(table=[0,68; 20,57; 40,50; 60,45; 80,43; 100,40; 120,38; 140,35; 160,32; 180,30; 200,28; 220,25], tableOnFile=false) "Polarization curve 5kW water cooled PEMFC (Zou et Kim)"
                                                                                                                                                                                                        annotation (Placement(transformation(
        extent={{-7,-7},{7,7}},
        rotation=180,
        origin={83,65})));
  Modelica.Blocks.Sources.CombiTimeTable PowerVoltageExp(table=[0,0; 10,-500; 20,-900; 30,-1400; 40,-1800; 50,-2200; 60,-2400; 70,-2700; 80,-3000; 90,-3400; 100,-3800; 110,-4000; 120,-4200; 130,-4500; 140,-4900; 150,-5000; 160,-5100; 170,-5200; 180,-5400; 190,-5500; 200,-5600; 210,-5500; 220,-5400], tableOnFile=false) "Polarization curve 5kW water cooled PEMFC (Zou et kim)" annotation (Placement(transformation(
        extent={{-7,-7},{7,7}},
        rotation=180,
        origin={83,43})));
  Modelica.Blocks.Sources.CombiTimeTable CurrentPolarization(table=[0,8; 0.1,8; 0.2,13; 0.3,18; 0.4,23; 0.5,28; 0.6,33; 0.7,38; 0.8,43; 0.9,48; 1.0,53; 1.1,58; 1.2,63; 1.3,68; 1.4,73; 1.5,78; 1.6,83; 1.7,88; 1.8,93; 1.9,98; 2.0,103; 2.1,108; 2.2,113; 2.3,118; 2.4,123; 2.5,128; 2.6,133; 2.7,138; 2.8,143; 2.9,148; 3.0,153; 3.1,158; 3.2,163; 3.3,168; 3.4,175], tableOnFile=false) "Laurencelle Ballard MK5-E polarization curve" annotation (Placement(transformation(
        extent={{-8,-8},{8,8}},
        rotation=0,
        origin={-82,30})));
  Modelica.Blocks.Sources.CombiTimeTable PolarizationV(
    table=[1,33.5; 2,32.8; 3,32.3; 4,31.8; 5,31.4; 6,30.9; 7,30.4; 8,29.9; 9,29.4; 10,28.9; 11,28.4; 12,27.9; 13,27.6; 14,27.3; 15,27; 16,26.7; 17,26.4; 18,26.1; 19,25.8; 20,25.5; 21,25.2; 22,24.9; 23,24.6; 24,24.3; 25,24],
    tableOnFile=false,
    timeScale=0.15) " Laurencelle Ballard MK5-E polarization curve Polarization curve 5kW air cooled PEMFC - voltage response" annotation (Placement(transformation(
        extent={{-7,-7},{7,7}},
        rotation=180,
        origin={83,19})));
  Modelica.Blocks.Sources.Constant CurrentSet(k=8) annotation (Placement(transformation(extent={{-6,74},{10,90}})));
  AirSupplySystem.AirCompressorSystem AirCompressorSystem annotation (Placement(transformation(extent={{24,-72},{44,-52}})));
  Modelica.Blocks.Sources.CombiTimeTable CurrentPolarization125A(table=[0,8; 0.1,8; 0.2,13; 0.3,18; 0.4,23; 0.5,28; 0.6,33; 0.7,38; 0.8,43; 0.9,48; 1.0,53; 1.1,58; 1.2,63; 1.3,68; 1.4,73; 1.5,78; 1.6,83; 1.7,88; 1.8,93; 1.9,98; 2.0,103; 2.1,108; 2.2,113; 2.3,118; 2.4,123], tableOnFile=false) "Laurencelle Ballard MK5-E polarization curve from 8 to 125A" annotation (Placement(transformation(
        extent={{-8,-8},{8,8}},
        rotation=0,
        origin={-34,28})));
  Modelica.Blocks.Sources.Ramp FullRangeCurrentRamp(
    height=222,
    duration=4.4,
    offset=8,
    startTime=0) "Laurencelle Ballard MK5-E polarization curve" annotation (Placement(transformation(extent={{-88,-48},{-72,-32}})));
  Modelica.Blocks.Sources.CombiTimeTable ExperimentalPolarizationCurve(
    table=[0,32.50; 1,31.7; 2,30.94; 3,30.25; 4,29.75; 5,29.46; 6,29.08; 7,28.7; 8,28.4; 9,28; 10,27.81; 11,27.55; 12,27.205; 13,27.04; 14,26.6; 15,26.42; 16,26.12; 17,25.8; 18,25.5; 19,25.1; 20,24.87; 21,24.51; 22,24.2; 23,23.7],
    tableOnFile=false,
    timeScale=0.1) " Laurencelle Ballard MK5-E polarization curve at 72°C for I = (8;125) A" annotation (Placement(transformation(
        extent={{-7,-7},{7,7}},
        rotation=180,
        origin={89,-9})));
equation

  connect(FC.feedh, SyngasSource.gas_a) annotation (Line(
      points={{-6,0},{-28,0},{-28,-1},{-32,-1}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.feeda, AirSource.gas_a) annotation (Line(
      points={{-6,-12},{-28,-12},{-28,-27},{-33,-27}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.drainh, SyngasSink.gas_a) annotation (Line(
      points={{14,0},{14,1},{54,1}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.draina, AirSink.gas_a) annotation (Line(
      points={{14,-12},{50,-12},{50,-28},{54,-28}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.lambda_H, lambdaController_PID.u1) annotation (Line(points={{9.2,-16},{9.2,-70},{-96,-70},{-96,-4},{-89.4,-4}}, color={0,0,127}));
  connect(lambdaController_PID.y, SyngasSource.m_flow) annotation (Line(points={{-70.6,0},{-54,0},{-54,3.8},{-48,3.8}},   color={0,0,127}));
  connect(FC.lambda_O, lambdaController_PID1.u1) annotation (Line(points={{-0.8,-16},{-0.8,-50},{-16.6,-50}}, color={0,0,127}));
  connect(lambdaController_PID1.y, AirSource.m_flow) annotation (Line(points={{-35.4,-54},{-52,-54},{-52,-32.4},{-46,-32.4}}, color={0,0,127}));
  connect(TempSet.y, coolingModel.T_environment) annotation (Line(points={{8.4,32},{20,32}}, color={0,0,127}));
  connect(FC.temperatureOut, coolingModel.T_op) annotation (Line(points={{1.8,-9},{1.8,-10},{-12,-10},{-12,40},{10,40},{10,38},{20,38}},   color={0,0,127}));
  connect(FC.heat, coolingModel.heatPortCooling) annotation (Line(points={{14.1,-9.1},{22,-9.1},{22,18},{16,18},{16,24.2},{20,24.2}}, color={191,0,0}));
  connect(lambdaController_PID1.y, AirCompressorSystem.AirMassFlowRateSetpoint) annotation (Line(points={{-35.4,-54},{-40,-54},{-40,-72},{18,-72},{18,-56},{24,-56}}, color={0,0,127}));
  connect(CurrentDensityStep.y, A_cell.u) annotation (Line(points={{-45.2,82},{-40,82},{-40,62},{-29,62},{-29,55}}, color={0,0,127}));
  connect(PolarizationCurrentRamp.y, FC.I_load) annotation (Line(points={{-49.2,30},{-46,30},{-46,12},{-22,12},{-22,-6.6},{-4.2,-6.6}}, color={0,0,127}));
  annotation (experiment(
      StopTime=4.2,
      Interval=1,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"));
end TestPEMFCexp;
