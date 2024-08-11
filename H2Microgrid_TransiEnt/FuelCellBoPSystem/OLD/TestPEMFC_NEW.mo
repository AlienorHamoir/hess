within H2Microgrid_TransiEnt.FuelCellBoPSystem.OLD;
model TestPEMFC_NEW "Example of a fuel cell in a domestic application that follows different power setpoints"

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
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{-92,84},{-78,96}})));

  parameter TransiEnt.Basics.Media.Gases.Gas_MoistAir Air=TransiEnt.Basics.Media.Gases.Gas_MoistAir() "Medium model of air" annotation (choicesAllMatching);

  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow AirSource(
    variable_m_flow=true,
    variable_xi=false,
    m_flow_const=0.001,
    T_const=23.5 + 273,
    medium=FC.Air)    annotation (Placement(transformation(
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
    medium=FC.Air)        annotation (Placement(transformation(
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
    xi_const={0,0,0,0,1,0})
                      annotation (Placement(transformation(
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
    xi_const={0,0,0,0,1,0})
                      annotation (Placement(transformation(extent={{-48,3},{-32,19}})));

  Modelica.Blocks.Sources.Ramp CurrentRamp(
    height=120,
    duration=4000,
    offset=1,
    startTime=100)
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
  AirSupplySystem.AirCompressorSystem airCompressorSystemModel annotation (Placement(transformation(extent={{-32,-94},{-12,-74}})));
  Modelica.Blocks.Sources.Ramp PowerRamp(
    height=4000,
    duration=1000,
    offset=500,
    startTime=10) annotation (Placement(transformation(extent={{-94,42},{-74,62}})));
  Modelica.Blocks.Sources.Constant PowerSet(k=600) annotation (Placement(transformation(extent={{4,66},{24,86}})));
  Modelica.Blocks.Sources.Step PowerStep(
    height=1000,
    offset=0,
    startTime=1000) annotation (Placement(transformation(extent={{-58,66},{-38,86}})));
  TransiEnt.Components.Electrical.FuelCellSystems.FuelCell.Controller.LambdaController lambdaController(m_flow_rampup=1e-8) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-70,-64})));
  Modelica.Thermal.FluidHeatFlow.Examples.Utilities.DoubleRamp Load1(
    startTime=1.25e3,
    interval=1.25e3,
    duration_1=1000,
    offset=0,
    height_1=5000,
    height_2=-4700,
    duration_2=1000)
                  annotation (Placement(transformation(extent={{62,48},{82,68}})));
  TransiEnt.Components.Electrical.FuelCellSystems.FuelCell.Controller.LambdaController lambdaController1(Lambda_H_target=2.05, m_flow_rampup=1e-8) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-26,-50})));
  Controller.PowerController powerController annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={8,32})));
  Controller.LambdaController_PID lambdaHController_PID(lambda_target=1.5, m_flow_rampup=1e-10) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={70,-74})));
equation

  connect(FC.feedh, SyngasSource.gas_a) annotation (Line(
      points={{-10,-2},{-26,-2},{-26,11},{-32,11}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.feeda, AirSource.gas_a) annotation (Line(
      points={{-10,-20},{-28,-20},{-28,-27},{-33,-27}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.drainh, SyngasSink.gas_a) annotation (Line(
      points={{22,-2},{42,-2},{42,9},{54,9}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.draina, AirSink.gas_a) annotation (Line(
      points={{22,-20},{38,-20},{38,-28},{54,-28}},
      color={255,170,85},
      thickness=0.5));

  connect(lambdaController.y, SyngasSource.m_flow) annotation (Line(points={{-80.8,-64},{-92,-64},{-92,15.8},{-48,15.8}}, color={0,0,127}));
  connect(FC.lambda_H, lambdaController.u1) annotation (Line(points={{14.32,-26},{14,-26},{14,-70},{-60,-70}}, color={0,0,127}));
  connect(FC.lambda_O, lambdaController1.u1) annotation (Line(points={{-1.68,-26},{-1.68,-56},{-16,-56}}, color={0,0,127}));
  connect(lambdaController1.y, AirSource.m_flow) annotation (Line(points={{-36.8,-50},{-52,-50},{-52,-32.4},{-46,-32.4}}, color={0,0,127}));
  connect(airCompressorSystemModel.AirMassFlowRateSetpoint, lambdaController1.y) annotation (Line(points={{-32,-78},{-46,-78},{-46,-50},{-36.8,-50}}, color={0,0,127}));
  connect(powerController.V_stack, FC.V_stack) annotation (Line(points={{17,37.4},{70,37.4},{70,-11},{22,-11}}, color={0,127,127}));
  connect(powerController.y, FC.I_load) annotation (Line(points={{-3,32},{-18,32},{-18,-11.9},{-7.12,-11.9}}, color={0,0,127}));
  connect(powerController.P, Load1.y) annotation (Line(points={{17,26},{90,26},{90,58},{83,58}}, color={0,127,127}));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
    experiment(
      StopTime=5000,
      Interval=1,
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
<p>Test environment for the modified PEMFC model</p>
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
<p>Based on TransiEnt library testing model structure. Compositions, medium, controllers and fuel cell models have been adapted.</p>
<h4><span style=\"color: #008000\">10. Version History</span></h4>
</html>"));
end TestPEMFC_NEW;
