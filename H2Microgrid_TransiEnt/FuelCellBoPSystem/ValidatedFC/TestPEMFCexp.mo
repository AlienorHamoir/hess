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
    offset=1,
    startTime=3000)
                  "To use as direct input to the fuel cell model, without power controller" annotation (Placement(transformation(extent={{34,66},{54,86}})));
  Modelica.Thermal.FluidHeatFlow.Examples.Utilities.DoubleRamp Load(
    startTime=1.25e4,
    interval=1.25e4,
    duration_1=1,
    offset=10,
    height_1=40,
    height_2=-15,
    duration_2=1)
                 annotation (Placement(transformation(extent={{-28,66},{-8,86}})));
  Modelica.Blocks.Sources.Ramp PowerRamp(
    height=30,
    duration=1000,
    offset=10,
    startTime=10) annotation (Placement(transformation(extent={{-94,42},{-74,62}})));
  Modelica.Blocks.Sources.Constant PowerSet(k=11)  annotation (Placement(transformation(extent={{4,66},{24,86}})));
  Modelica.Blocks.Sources.Step PowerStep(
    height=50,
    offset=0,
    startTime=500)  annotation (Placement(transformation(extent={{-58,66},{-38,86}})));
  PEMFC_KhanIqbal FC(usePowerPort=false, useHeatPort=false)
                     annotation (Placement(transformation(extent={{-6,-18},{14,2}})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{-90,82},{-74,98}})));
  TransiEnt.Components.Electrical.FuelCellSystems.FuelCell.Controller.LambdaController lambdaController(m_flow_rampup=1e-8) annotation (Placement(transformation(extent={{-90,4},{-70,24}})));
  TransiEnt.Components.Electrical.FuelCellSystems.FuelCell.Controller.LambdaController lambdaControllerO2(
    Lambda_H_target=2.05,
    m_flow_rampup=1e-4,
    k=1e-2) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-18,-44})));
equation

  connect(FC.feedh, SyngasSource.gas_a) annotation (Line(
      points={{-6,-2},{-28,-2},{-28,11},{-32,11}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.feeda, AirSource.gas_a) annotation (Line(
      points={{-6,-14},{-28,-14},{-28,-27},{-33,-27}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.drainh, SyngasSink.gas_a) annotation (Line(
      points={{14,-2},{14,9},{54,9}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.draina, AirSink.gas_a) annotation (Line(
      points={{14,-14},{50,-14},{50,-28},{54,-28}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.lambda_H, lambdaController.u1) annotation (Line(points={{9.2,-18},{9.2,-26},{8,-26},{8,-96},{-94,-96},{-94,20},{-90,20}}, color={0,0,127}));
  connect(lambdaController.y, SyngasSource.m_flow) annotation (Line(points={{-69.2,14},{-52,14},{-52,15.8},{-48,15.8}}, color={0,0,127}));
  connect(FC.lambda_O, lambdaControllerO2.u1) annotation (Line(points={{-0.8,-18},{-0.8,-50},{-8,-50}}, color={0,0,127}));
  connect(lambdaControllerO2.y, AirSource.m_flow) annotation (Line(points={{-28.8,-44},{-56,-44},{-56,-32.4},{-46,-32.4}}, color={0,0,127}));
  connect(Load.y, FC.I_load) annotation (Line(points={{-7,76},{-4,76},{-4,6},{-14,6},{-14,-8.6},{-4.2,-8.6}}, color={0,0,127}));
  annotation (experiment(
      StartTime=1,
      StopTime=40000,
      Interval=1,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"));
end TestPEMFCexp;
