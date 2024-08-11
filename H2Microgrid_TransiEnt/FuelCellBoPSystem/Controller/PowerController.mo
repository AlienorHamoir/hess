within H2Microgrid_TransiEnt.FuelCellBoPSystem.Controller;
model PowerController "Controller for power and current output in Fuel Cell applications"

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

  // _____________________________________________
  //
  //          Imports and Class Hierarchy
  // _____________________________________________

  extends TransiEnt.Basics.Icons.Controller;

  // _____________________________________________
  //
  //             Visible Parameters
  // _____________________________________________

  parameter Modelica.Units.SI.Current i_max=150 "Maximum current value";
  parameter Real k = 1 "Proportional controller gain";

  // _____________________________________________
  //
  //                  Interfaces
  // _____________________________________________

  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn P "Input for power setpoint" annotation (Placement(
        transformation(
        rotation=180,
        extent={{10,-10},{-10,10}},
        origin={-100,74}), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-90,60})));
  Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(transformation(
          rotation=180,
                      extent={{10,-10},{-10,10}},
        origin={104,0}), iconTransformation(extent={{100,-10},{120,10}})));
  TransiEnt.Basics.Interfaces.Electrical.VoltageIn V_stack "Input for stack voltage" annotation (Placement(
        transformation(
        rotation=0,
        extent={{-10,-10},{10,10}},
        origin={-100,-66}), iconTransformation(
        extent={{-10,10},{10,-10}},
        rotation=0,
        origin={-90,-54})));

  // _____________________________________________
  //
  //           Instances of other Classes
  // _____________________________________________

  Modelica.Blocks.Math.Gain Gain(k=k) annotation (Placement(transformation(
        extent={{6,-6.5},{-6,6.5}},
        rotation=180,
        origin={-42,24.5})));

  Modelica.Blocks.Math.Division PowerbyVoltage_divison annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=180,
        origin={5,0})));

  Modelica.Blocks.Nonlinear.Limiter constantSaturation(uMax=i_max, uMin=0.1)
    annotation (Placement(transformation(extent={{48,-12},{72,12}})));

  Modelica.Blocks.Nonlinear.Limiter preventDivisionByZero(uMax=150,  uMin=0.1)
    annotation (Placement(transformation(extent={{-50,-38},{-32,-20}})));

equation
  // _____________________________________________
  //
  //               Connect Statements
  // _____________________________________________

  connect(constantSaturation.y, y) annotation (Line(
      points={{73.2,0},{104,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(Gain.y, PowerbyVoltage_divison.u1) annotation (Line(
      points={{-35.4,24.5},{-32,24.5},{-32,6},{-7,6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(preventDivisionByZero.y, PowerbyVoltage_divison.u2) annotation (Line(
      points={{-31.1,-29},{-31.1,-6},{-7,-6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(V_stack, preventDivisionByZero.u) annotation (Line(points={{-100,-66},{-64,-66},{-64,-29},{-51.8,-29}},
                                                                                              color={0,127,127}));
  connect(P, Gain.u) annotation (Line(points={{-100,74},{-64,74},{-64,24.5},{-49.2,24.5}},   color={0,127,127}));
  connect(PowerbyVoltage_divison.y, constantSaturation.u) annotation (Line(points={{16,0},{45.6,0}}, color={0,0,127}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})),           Icon(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={Text(
          extent={{-126,94},{-76,74}},
          lineColor={0,0,255},
          textString="deltaQ"), Text(
          extent={{-128,-20},{-78,-40}},
          lineColor={0,0,255},
          textString="Vcell")}),
                Documentation(info="<html>
<h4><span style=\"color: #008000\">1. Purpose of model</span></h4>
<p>Controller for power output in Fuel Cell applications</p>
<p>Based on input power setpoint and stack voltage, the current setpoint for the fuel cell is calculated. </p>
<h4><span style=\"color: #008000\">2. Level of detail, physical effects considered, and physical insight</span></h4>
<p>It has been adapted from the initial model by differentiating several cases:</p>
<ol>
<li>If stack voltage is larger than 0, I = P / V</li>
<li>If power setpoint is smaller than minimum allowed power and stack voltage is smaller or equal to 0, the current is 0</li>
<li>If power setpoint is larger than minimum allowed power and stack voltage is smaller or equal to 0, the current equals rampup current</li>
</ol>
<h4><span style=\"color: #008000\">3. Limits of validity </span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">4. Interfaces</span></h4>
<p>Modelica RealInput: electric voltage in V</p>
<p>Modelica RealInput: electric power difference in W</p>
<p>Modelica RealOutput: y</p>
<h4><span style=\"color: #008000\">5. Nomenclature</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">6. Governing Equations</span></h4>
<p>I = P / V</p>
<h4><span style=\"color: #008000\">7. Remarks for Usage</span></h4>
<p>Adapt current, voltage and power operating ranges and characteristics</p>
<h4><span style=\"color: #008000\">8. Validation</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">9. References</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">10. Version History</span></h4>
<p>Model created by Simon Weilbach (simon.weilbach@tuhh.de) on 01.10.2014</p>
<p>Model revised by Pascal Dubucq (dubucq@tuhh.de) on 01.10.2015</p>
<p>Quality check (Code conventions) by Rebekka Denninger on 01.10.2016</p>
<p>Model adapted by Ali&eacute;nor Hamoir in June 2024</p>
</html>"));
end PowerController;
