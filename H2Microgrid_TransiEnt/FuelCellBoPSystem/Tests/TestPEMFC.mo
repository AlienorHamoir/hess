within H2Microgrid_TransiEnt.FuelCellBoPSystem.Tests;
model TestPEMFC "Example of a fuel cell in a domestic application that follows load such that power grid consumption is minimized"

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
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{-90,80},{-70,100}})));

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

  FuelCell.PEMFC FC(
    Syngas=TransiEnt.Basics.Media.Gases.Gas_VDIWA_SG7_var(),
    I_shutdown=10,
    T_nom(displayUnit="K") = 350,
    T_stack_max(displayUnit="K") = 355,
    T_cool_set(displayUnit="K") = 350,
    usePowerPort=false,
    useHeatPort=true,
    I(start=10))                  annotation (Placement(transformation(extent={{-10,-26},{22,4}})));
  Modelica.Blocks.Sources.Ramp ramp(
    height=90,
    duration=400,
    offset=10,
    startTime=10) annotation (Placement(transformation(extent={{-16,58},{4,78}})));
  FuelCell.Controller.LambdaController_PID lambdaHController_PID(lambda_target=1.5, m_flow_rampup=1e-6) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-72,-66})));
  FuelCell.Controller.LambdaController_PID lambdaOController_PID(lambda_target=2.05, m_flow_rampup=2e-6)
                                                                 "Controller that outputs the required air mass flow rate to meet OER (oxygen excess ratio) target " annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-18,-50})));
  CoolingSystem.HeatPortCooling.CoolingModel coolingModel annotation (Placement(transformation(extent={{70,48},{90,68}})));
  Modelica.Thermal.FluidHeatFlow.Examples.Utilities.DoubleRamp Load(
    startTime=1.25e4,
    interval=1.25e4,
    duration_1=1,
    offset=10,
    height_1=40,
    height_2=-15,
    duration_2=1)
                 annotation (Placement(transformation(extent={{-46,68},{-26,88}})));
  AirSupplySystem.AirCompressorSystemModel airCompressorSystemModel annotation (Placement(transformation(extent={{-26,-90},{-6,-70}})));
  Modelica.Blocks.Sources.Ramp ramp1(
    height=4000,
    duration=1000,
    offset=500,
    startTime=10) annotation (Placement(transformation(extent={{-94,42},{-74,62}})));
  Modelica.Blocks.Sources.Constant const(k=600) annotation (Placement(transformation(extent={{28,68},{48,88}})));
  FuelCell.Controller.PowerConverter powerConverter(i_max=150) annotation (Placement(transformation(extent={{-52,26},{-32,46}})));
  Modelica.Blocks.Sources.Step step(
    height=1000,
    offset=0,
    startTime=1000) annotation (Placement(transformation(extent={{-118,14},{-98,34}})));
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

  connect(FC.lambda_H, lambdaHController_PID.u1) annotation (Line(points={{14.32,-26},{14.32,-62.4},{-61.8,-62.4}},
                                                                                                              color={0,0,127}));
  connect(lambdaHController_PID.y, SyngasSource.m_flow) annotation (Line(points={{-82.8,-66},{-88,-66},{-88,16},{-58,16},{-58,15.8},{-48,15.8}}, color={0,0,127}));
  connect(FC.lambda_O, lambdaOController_PID.u1) annotation (Line(points={{-1.68,-26},{-1.68,-36},{-2,-36},{-2,-46},{-6,-46},{-6,-46.4},{-7.8,-46.4}},color={0,0,127}));
  connect(FC.heat, coolingModel.heatPortCooling) annotation (Line(points={{22.16,-15.65},{70,-15.65},{70,48.8}},         color={191,0,0}));
  connect(FC.temperatureOut, coolingModel.T_op) annotation (Line(points={{-1.04,-15.5},{10,-15.5},{10,36},{64,36},{64,64.6},{68.8,64.6}},   color={0,0,127}));
  connect(lambdaOController_PID.y, AirSource.m_flow) annotation (Line(points={{-28.8,-50},{-52,-50},{-52,-32.4},{-46,-32.4}}, color={0,0,127}));
  connect(lambdaOController_PID.y, airCompressorSystemModel.AirMassFlowRateSetpoint) annotation (Line(points={{-28.8,-50},{-34,-50},{-34,-73.6},{-26.8,-73.6}}, color={0,0,127}));
  connect(FC.V_stack, powerConverter.V_stack) annotation (Line(
      points={{22,-11},{22,-12},{30,-12},{30,50},{-58,50},{-58,30.6},{-51,30.6}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(powerConverter.y, FC.I_load) annotation (Line(points={{-31,36},{-22,36},{-22,-11.9},{-7.12,-11.9}}, color={0,0,127}));
  connect(ramp1.y, powerConverter.P) annotation (Line(points={{-73,52},{-56,52},{-56,42},{-51,42}}, color={0,0,127}));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
    experiment(
      StopTime=2000,
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
<p>Test environment for the PEM model</p>
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
end TestPEMFC;
