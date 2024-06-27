within H2Microgrid_TransiEnt.FuelCellBoPSystem.AirSupplySystem;
model TestAirCompressor "Example of a fuel cell in a domestic application that follows load such that power grid consumption is minimized"

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
    duration_1=50,
    duration_2=55,
    offset=1000,
    height_1=200,
    height_2=-400)
                 annotation (Placement(transformation(extent={{-42,70},{-22,90}})));
//     xi_const={0.01,0.7},
        //     xi_const={1,0},

  Modelica.Blocks.Sources.Ramp ramp(
    height=0.4000,
    duration=500,
    offset=0,
    startTime=10) annotation (Placement(transformation(extent={{-62,8},{-42,28}})));
  AirSupplySystem.AirCompressor airCompressor annotation (Placement(transformation(extent={{6,6},{26,26}})));
  Modelica.Blocks.Sources.Step step(height=0.005, startTime=15) annotation (Placement(transformation(extent={{-54,-32},{-34,-12}})));
equation

  connect(step.y, airCompressor.AirMassFlowRateSetpoint) annotation (Line(points={{-33,-22},{-2,-22},{-2,24},{5,24}}, color={0,0,127}));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
    experiment(
      StopTime=1000,
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
end TestAirCompressor;
