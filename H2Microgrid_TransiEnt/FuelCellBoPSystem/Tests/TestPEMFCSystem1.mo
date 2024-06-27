within H2Microgrid_TransiEnt.FuelCellBoPSystem.Tests;
model TestPEMFCSystem1 "Example of a fuel cell in a domestic application that follows load such that power grid consumption is minimized"

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

  parameter H2Microgrid_TransiEnt.FuelCellBoPSystem.FuelCell.Physics.Gas_VDIWA_H2_var Syngas=H2Microgrid_TransiEnt.FuelCellBoPSystem.FuelCell.Physics.Gas_VDIWA_H2_var() "Medium model H2" annotation (Dialog(group="Fundamental Definitions"));

  parameter TransiEnt.Basics.Media.Gases.Gas_MoistAir Air=TransiEnt.Basics.Media.Gases.Gas_MoistAir() "Medium model of air" annotation (choicesAllMatching);

  Modelica.Thermal.FluidHeatFlow.Examples.Utilities.DoubleRamp Load(
    startTime=5,
    interval=5,
    duration_1=5,
    duration_2=5,
    offset=1000,
    height_1=200,
    height_2=-400)
                 annotation (Placement(transformation(extent={{-30,56},{-10,76}})));
TransiEnt.Components.Electrical.FuelCellSystems.FuelCell.Controller.PowerController PowerController(k=1)    annotation (Placement(transformation(rotation=0, extent={{-34,18},{-54,38}})));
//     xi_const={0.01,0.7},
  ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi AirSink(
    variable_p=false,
    variable_T=true,
    variable_xi=false,
    p_const=1e5,
    T_const=80 + 273.15,
    medium=FC.Air)        annotation (Placement(transformation(
        extent={{-7,-8},{7,8}},
        rotation=180,
        origin={59,-50})));
        //     xi_const={1,0},
  ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi SyngasSink(
    variable_p=false,
    variable_T=true,
    variable_xi=false,
    p_const=1e5,
    medium=FC.Syngas,
    T_const=80 + 273.15,
    xi_const={0,0,0.003,0.995,0.001,0})
                      annotation (Placement(transformation(
        extent={{-6,-9},{6,9}},
        rotation=180,
        origin={58,-13})));

  FuelCell.PEMFC FC(
    Syngas=TransiEnt.Basics.Media.Gases.Gas_VDIWA_SG7_var(),
    m=1,
    cp=35000,
    no_Cells=35,
    T_nom=308.15,
    T_stack_max=313.15,
    T_cool_set=303.15,
    T_stack(start=298.15),
    usePowerPort=false,
    useHeatPort=true)             annotation (Placement(transformation(extent={{-16,-48},{16,-18}})));
  Modelica.Blocks.Sources.Ramp ramp(
    height=4000,
    duration=500,
    offset=500,
    startTime=10) annotation (Placement(transformation(extent={{14,54},{34,74}})));
  FuelCell.Controller.LambdaController_PID lambdaHController_PID(lambda_target=1.5, m_flow_rampup=1e-6) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-72,-88})));
  FuelCell.Controller.LambdaController_PID lambdaOController_PID annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-26,-72})));
  CoolingSystem.HeatPortCooling.CoolingModel coolingModel annotation (Placement(transformation(extent={{68,-2},{88,18}})));
  TransiEnt.Basics.Adapters.Gas.Real_ph_to_Ideal_pT real_ph_to_Ideal_pT(redeclare TransiEnt.Basics.Media.Gases.VLE_VDIWA_H2 real, redeclare TransiEnt.Basics.Media.Gases.Gas_VDIWA_SG7_var ideal) annotation (Placement(transformation(extent={{-48,-30},{-34,-16}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_Txim_flow boundary_Txim_flow(
    medium=TransiEnt.Basics.Media.Gases.VLE_VDIWA_H2(),
    variable_m_flow=true,
    T_const=313.15) annotation (Placement(transformation(extent={{-80,-30},{-66,-16}})));
  AirSupplySystem.AirCompressorSystemModel airCompressorSystemModel annotation (Placement(transformation(
        extent={{-15,-9},{15,9}},
        rotation=0,
        origin={47,-77})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow AirSource(
    variable_m_flow=true,
    variable_xi=false,
    m_flow_const=0.001,
    T_const=25 + 273,
    medium=FC.Air)    annotation (Placement(transformation(
        extent={{6.5,-9},{-6.5,9}},
        rotation=180,
        origin={-51.5,-47})));
equation

  connect(FC.drainh, SyngasSink.gas_a) annotation (Line(
      points={{16,-24},{40,-24},{40,-13},{52,-13}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.draina, AirSink.gas_a) annotation (Line(
      points={{16,-42},{36,-42},{36,-50},{52,-50}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.feeda, AirSource.gas_a) annotation (Line(
      points={{-16,-42},{-34,-42},{-34,-47},{-45,-47}},
      color={255,170,85},
      thickness=0.5));
  connect(PowerController.V_cell, FC.V_stack) annotation (Line(points={{-35,22.6},{26,22.6},{26,-33},{16,-33}},
                                                                                                    color={0,0,127}));
  connect(PowerController.y, FC.I_load) annotation (Line(points={{-55,28},{-60,28},{-60,-33.9},{-13.12,-33.9}},            color={0,0,127}));

  connect(ramp.y, PowerController.deltaP) annotation (Line(points={{35,64},{40,64},{40,34},{-35,34}}, color={0,0,127}));
  connect(FC.lambda_H, lambdaHController_PID.u1) annotation (Line(points={{8.32,-48},{8.32,-84.4},{-61.8,-84.4}},
                                                                                                              color={0,0,127}));
  connect(FC.lambda_O, lambdaOController_PID.u1) annotation (Line(points={{-7.68,-48},{-7.68,-58},{-8,-58},{-8,-68},{-12,-68},{-12,-68.4},{-15.8,-68.4}},
                                                                                                                                                      color={0,0,127}));
  connect(FC.heat, coolingModel.heatPortCooling) annotation (Line(points={{16,-36.9},{68,-36.9},{68,-1.2}},              color={191,0,0}));
  connect(FC.temperatureOut, coolingModel.T_op) annotation (Line(points={{-7.04,-37.5},{4,-37.5},{4,14},{66,14},{66,14.6},{66.8,14.6}},     color={0,0,127}));
  connect(FC.temperatureOut, AirSink.T) annotation (Line(points={{-7.04,-37.5},{4,-37.5},{4,14},{58,14},{58,22},{92,22},{92,-50},{66,-50}}, color={0,0,127}));
  connect(FC.temperatureOut, SyngasSink.T) annotation (Line(points={{-7.04,-37.5},{4,-37.5},{4,14},{58,14},{58,22},{92,22},{92,-13},{64,-13}}, color={0,0,127}));
  connect(real_ph_to_Ideal_pT.gasPortOut, FC.feedh) annotation (Line(
      points={{-34,-23},{-24,-23},{-24,-24},{-16,-24}},
      color={255,170,85},
      thickness=1.5));
  connect(real_ph_to_Ideal_pT.gasPortIn, boundary_Txim_flow.gasPort) annotation (Line(
      points={{-48,-23},{-66,-23}},
      color={255,255,0},
      thickness=1.5));
  connect(boundary_Txim_flow.m_flow, lambdaHController_PID.y) annotation (Line(points={{-81.4,-18.8},{-92,-18.8},{-92,-90},{-82.8,-90},{-82.8,-88}}, color={0,0,127}));
  connect(lambdaOController_PID.y, airCompressorSystemModel.massFlowRateIn) annotation (Line(points={{-36.8,-72},{-40,-72},{-40,-86},{26,-86},{26,-70},{41.6,-70},{41.6,-67.82}}, color={0,0,127}));
  connect(lambdaOController_PID.y, AirSource.m_flow) annotation (Line(points={{-36.8,-72},{-64,-72},{-64,-52.4},{-58,-52.4}}, color={0,0,127}));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
    experiment(
      StopTime=6000,
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
end TestPEMFCSystem1;
