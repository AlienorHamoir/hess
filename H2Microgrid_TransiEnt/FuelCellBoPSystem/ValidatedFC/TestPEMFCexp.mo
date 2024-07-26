within H2Microgrid_TransiEnt.FuelCellBoPSystem.ValidatedFC;
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
        origin={60,9})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow SyngasSource(
    variable_T=false,
    m_flow_const=5.1e-2,
    variable_m_flow=true,
    variable_xi=false,
    T_const=40 + 273.15,
    medium=FC.Syngas,
    xi_const={0,0,0,0,1,0}) annotation (Placement(transformation(extent={{-48,3},{-32,19}})));

  Modelica.Blocks.Sources.Ramp CurrentRamp(
    height=20,
    duration=1,
    offset=0.1,
    startTime=1000)
                  "To use as direct input to the fuel cell model, without power controller" annotation (Placement(transformation(extent={{34,66},{54,86}})));
  Modelica.Thermal.FluidHeatFlow.Examples.Utilities.DoubleRamp DynamicCurrentDensity(
    startTime=13,
    interval=26,
    duration_1=1,
    offset=0.015,
    height_1=0.095,
    height_2=-0.095,
    duration_2=1) "Test current density for dynamical model validation" annotation (Placement(transformation(extent={{-28,66},{-8,86}})));
  Modelica.Blocks.Sources.Ramp PowerRamp(
    height=100,
    duration=1,
    offset=10,
    startTime=1000)
                  annotation (Placement(transformation(extent={{-88,72},{-68,92}})));
  Modelica.Blocks.Sources.Constant PowerSet(k=15 + 273.15)
                                                   annotation (Placement(transformation(extent={{4,66},{24,86}})));
  Modelica.Blocks.Sources.Step CurrentDensityStep(
    height=-0.49,
    offset=0.5,
    startTime=100)
                  annotation (Placement(transformation(extent={{-58,66},{-38,86}})));
  PEMFC_KhanIqbal FC(
    p_Anode=200000,     usePowerPort=false,
    useHeatPort=true,
    T_stack(start=343.15))
                     annotation (Placement(transformation(extent={{-6,-16},{14,4}})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{72,-94},{88,-78}})));
  Modelica.Blocks.Math.Gain A_cell(k=232) annotation (Placement(transformation(
        extent={{-5,-5},{5,5}},
        rotation=-90,
        origin={-19,45})));
  Modelica.Blocks.Sources.CombiTimeTable CellVoltageExp(table=[0,0.995; 13,0.995; 14,0.82; 16,0.8; 28,0.8; 40,0.8; 41,0.97; 42,0.993; 43,0.995; 70,0.995], tableOnFile=false) "cell dynamic voltage response (Kan et Iqbal)" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={76,72})));
  FuelCell.Controller.LambdaController_PID lambdaController_PID(lambda_target=1.5, m_flow_rampup=1e-8) annotation (Placement(transformation(extent={{-90,-10},{-70,10}})));
  FuelCell.Controller.LambdaController_PID lambdaController_PID1(lambda_target=2.05, m_flow_rampup=1e-8) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-26,-54})));
  Modelica.Thermal.FluidHeatFlow.Examples.Utilities.DoubleRamp LoadChange(
    startTime=10,
    interval=10,
    duration_1=0.1,
    offset=47.55,
    height_1=-37.55,
    height_2=75.15,
    duration_2=0.1) "Test current density for dynamical model validation" annotation (Placement(transformation(extent={{-92,34},{-72,54}})));
  CoolingSystem.HeatPortCooling.CoolingModel coolingModel(k_p=1000, tau_i=0.01) annotation (Placement(transformation(extent={{20,22},{40,42}})));
equation

  connect(FC.feedh, SyngasSource.gas_a) annotation (Line(
      points={{-6,0},{-28,0},{-28,11},{-32,11}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.feeda, AirSource.gas_a) annotation (Line(
      points={{-6,-12},{-28,-12},{-28,-27},{-33,-27}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.drainh, SyngasSink.gas_a) annotation (Line(
      points={{14,0},{14,9},{54,9}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.draina, AirSink.gas_a) annotation (Line(
      points={{14,-12},{50,-12},{50,-28},{54,-28}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.lambda_H, lambdaController_PID.u1) annotation (Line(points={{9.2,-16},{9.2,-70},{-96,-70},{-96,-4},{-89.4,-4}}, color={0,0,127}));
  connect(lambdaController_PID.y, SyngasSource.m_flow) annotation (Line(points={{-70.6,0},{-54,0},{-54,15.8},{-48,15.8}}, color={0,0,127}));
  connect(FC.lambda_O, lambdaController_PID1.u1) annotation (Line(points={{-0.8,-16},{-0.8,-50},{-16.6,-50}}, color={0,0,127}));
  connect(lambdaController_PID1.y, AirSource.m_flow) annotation (Line(points={{-35.4,-54},{-52,-54},{-52,-32.4},{-46,-32.4}}, color={0,0,127}));
  connect(PowerSet.y, coolingModel.T_environment) annotation (Line(points={{25,76},{30,76},{30,46},{12,46},{12,32},{20,32}}, color={0,0,127}));
  connect(FC.temperatureOut, coolingModel.T_op) annotation (Line(points={{-0.4,-9},{-0.4,-10},{-12,-10},{-12,36},{10,36},{10,38},{20,38}}, color={0,0,127}));
  connect(FC.heat, coolingModel.heatPortCooling) annotation (Line(points={{14.1,-9.1},{22,-9.1},{22,18},{16,18},{16,24.2},{20,24.2}}, color={191,0,0}));
  connect(A_cell.y, FC.I_load) annotation (Line(points={{-19,39.5},{-20,39.5},{-20,-6.6},{-4.2,-6.6}}, color={0,0,127}));
  connect(CurrentDensityStep.y, A_cell.u) annotation (Line(points={{-37,76},{-32,76},{-32,56},{-19,56},{-19,51}}, color={0,0,127}));
  annotation (experiment(
      StartTime=1,
      StopTime=4000,
      Interval=1,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"));
end TestPEMFCexp;
