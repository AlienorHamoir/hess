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
                 annotation (Placement(transformation(extent={{-236,44},{-216,64}})));
TransiEnt.Components.Electrical.FuelCellSystems.FuelCell.Controller.PowerController PowerController(k=1, i_min=20)
                                                                                                            annotation (Placement(transformation(rotation=0, extent={{-34,40},{-54,60}})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow AirSource(
    variable_m_flow=true,
    variable_xi=false,
    m_flow_const=0.001,
    T_const=25 + 273,
    medium=FC.Air)    annotation (Placement(transformation(
        extent={{6.5,-9},{-6.5,9}},
        rotation=180,
        origin={-39.5,-27})));
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
        origin={61,-28})));
        //     xi_const={1,0},
  ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi SyngasSink(
    variable_p=false,
    variable_T=true,
    variable_xi=false,
    p_const=1e5,
    medium=FC.Syngas,
    T_const=80 + 273.15,
    xi_const={0,0,0,0.005,0.99,0})
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
    xi_const={0,0,0,0.005,0.99,0})
                      annotation (Placement(transformation(extent={{-48,3},{-32,19}})));

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
    useHeatPort=true)             annotation (Placement(transformation(extent={{-10,-26},{22,4}})));
  Modelica.Blocks.Sources.Ramp ramp(
    height=4000,
    duration=500,
    offset=500,
    startTime=10) annotation (Placement(transformation(extent={{12,72},{32,92}})));
  FuelCell.Controller.LambdaController_PID lambdaHController_PID(lambda_target=1.5, m_flow_rampup=1e-6) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-72,-66})));
  FuelCell.Controller.LambdaController_PID lambdaOController_PID annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-18,-50})));
  CoolingSystem.HeatPortCooling.CoolingModel coolingModel annotation (Placement(transformation(extent={{70,48},{90,68}})));
  AirSupplySystem.AirCompressor airCompressor annotation (Placement(transformation(
        extent={{-10,-9},{10,9}},
        rotation=-90,
        origin={-39,-84})));
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
  connect(PowerController.V_cell, FC.V_stack) annotation (Line(points={{-35,44.6},{32,44.6},{32,-11},{22,-11}},
                                                                                                    color={0,0,127}));
  connect(PowerController.y, FC.I_load) annotation (Line(points={{-55,50},{-54,50},{-54,-11.9},{-7.12,-11.9}},             color={0,0,127}));

  connect(ramp.y, PowerController.deltaP) annotation (Line(points={{33,82},{38,82},{38,56},{-35,56}}, color={0,0,127}));
  connect(FC.lambda_H, lambdaHController_PID.u1) annotation (Line(points={{14.32,-26},{14.32,-62.4},{-61.8,-62.4}},
                                                                                                              color={0,0,127}));
  connect(lambdaHController_PID.y, SyngasSource.m_flow) annotation (Line(points={{-82.8,-66},{-88,-66},{-88,16},{-58,16},{-58,15.8},{-48,15.8}}, color={0,0,127}));
  connect(FC.lambda_O, lambdaOController_PID.u1) annotation (Line(points={{-1.68,-26},{-1.68,-36},{-2,-36},{-2,-46},{-6,-46},{-6,-46.4},{-7.8,-46.4}},color={0,0,127}));
  connect(FC.heat, coolingModel.heatPortCooling) annotation (Line(points={{22.16,-15.65},{70,-15.65},{70,48.8}},         color={191,0,0}));
  connect(FC.temperatureOut, coolingModel.T_op) annotation (Line(points={{-1.04,-15.5},{10,-15.5},{10,36},{64,36},{64,64.6},{68.8,64.6}},   color={0,0,127}));
  connect(FC.temperatureOut, AirSink.T) annotation (Line(points={{-1.04,-15.5},{10,-15.5},{10,36},{98,36},{98,-28},{68,-28}}, color={0,0,127}));
  connect(FC.temperatureOut, SyngasSink.T) annotation (Line(points={{-1.04,-15.5},{10,-15.5},{10,36},{98,36},{98,9},{66,9}}, color={0,0,127}));
  connect(lambdaOController_PID.y, airCompressor.AirMassFlowRateSetpoint) annotation (Line(points={{-28.8,-50},{-32,-50},{-32,-66},{-31.8,-66},{-31.8,-73}}, color={0,0,127}));
  connect(lambdaOController_PID.y, AirSource.m_flow) annotation (Line(points={{-28.8,-50},{-52,-50},{-52,-32.4},{-46,-32.4}}, color={0,0,127}));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
    experiment(
      StopTime=1500,
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
end TestPEMFC;
