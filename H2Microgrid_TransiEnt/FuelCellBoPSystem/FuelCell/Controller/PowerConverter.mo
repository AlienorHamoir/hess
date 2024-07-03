within H2Microgrid_TransiEnt.FuelCellBoPSystem.FuelCell.Controller;
model PowerConverter "Controller for power and current output in Fuel Cell applications"

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

  parameter Modelica.Units.SI.Current i_max=90 "Maximum current value";
  parameter Modelica.Units.SI.Current i_rampup=10 "Rampup current value";
  parameter Real k = 1 "Proportional controller gain";
  parameter Modelica.Units.SI.Power P_min=367.5 "Lower power Limit"; // computed from I_shutdown = 10 A and V_shutdown = 36.75 V

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
        origin={110,0}), iconTransformation(extent={{100,-10},{120,10}})));
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
        origin={-40,8.5})));

  Modelica.Blocks.Math.Division PowerbyVoltage_divison annotation (Placement(transformation(
        extent={{13.5,13},{-13.5,-13}},
        rotation=180,
        origin={0.5,-1})));

  Modelica.Blocks.Nonlinear.Limiter constantSaturation(uMax=i_max, uMin=0)
    annotation (Placement(transformation(extent={{68,-12},{92,12}})));

  Modelica.Blocks.Nonlinear.Limiter preventDivisionByZero(uMax=150,  uMin=0.1)
    annotation (Placement(transformation(extent={{-56,-54},{-38,-36}})));

equation
  // _____________________________________________
  //
  //               Connect Statements
  // _____________________________________________

  connect(constantSaturation.y, y) annotation (Line(
      points={{93.2,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(Gain.y, PowerbyVoltage_divison.u1) annotation (Line(
      points={{-33.4,8.5},{-26,8.5},{-26,6.8},{-15.7,6.8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(preventDivisionByZero.y, PowerbyVoltage_divison.u2) annotation (Line(
      points={{-37.1,-45},{-37.1,-46},{-26,-46},{-26,-8.8},{-15.7,-8.8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(V_stack, preventDivisionByZero.u) annotation (Line(points={{-100,-66},{-62,-66},{-62,-45},{-57.8,-45}},
                                                                                              color={0,127,127}));
  connect(P, Gain.u) annotation (Line(points={{-100,74},{-60,74},{-60,8.5},{-47.2,8.5}},     color={0,127,127}));
  connect(PowerbyVoltage_divison.y, constantSaturation.u) annotation (Line(points={{15.35,-1},{15.35,0},{65.6,0}}, color={0,0,127}));
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
<p>Basic model that divides input power by input voltage, to obtain a current setpoint. </p>
<h4><span style=\"color: #008000\">2. Level of detail, physical effects considered, and physical insight</span></h4>
<p>(no remarks)</p>
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
<p>Current operating range must be adapted for the application. It works if a superior controller only delivers power setpoint </p>
<h4><span style=\"color: #008000\">8. Validation</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">9. References</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">10. Version History</span></h4>
<p>Model created by Ali&eacute;nor Hamoir in June 2024</p>
</html>"));
end PowerConverter;
