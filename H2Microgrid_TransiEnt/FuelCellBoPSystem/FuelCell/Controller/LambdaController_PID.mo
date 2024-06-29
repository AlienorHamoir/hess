within H2Microgrid_TransiEnt.FuelCellBoPSystem.FuelCell.Controller;
model LambdaController_PID "PID Controller for Lambda in Fuel Cell Applications - it outputs the required air or hydrogen mass flow rate to meet OER/HER (oxygen/hydrogen excess ratio) target"

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

  parameter Real lambda_target=2.5 " excess ratio target - value for O2";
  parameter Modelica.Units.SI.MassFlowRate m_flow_rampup=2e-4 "Minimum mass flow";

  Boolean lambda_control;

     //Oxygen mass flow rate PID Controller Parameters
  parameter Modelica.Units.SI.Time tau_i=0.1 "1/tau_i for system PID integrator gain ";
  parameter Real k_p=0.001 "gain, cooling system PID proportional control ";
  parameter Real N_i=0.5 "gain of anti-windup compensation ";

  // _____________________________________________
  //
  //                  Interfaces
  // _____________________________________________

  Modelica.Blocks.Interfaces.RealInput u1 annotation (Placement(transformation(
          rotation=0, extent={{-112,-46},{-92,-26}})));
  Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(transformation(
          rotation=0, extent={{98,-10},{118,10}})));

  // _____________________________________________
  //
  //           Instances of other Classes
  // _____________________________________________

  Modelica.Blocks.Sources.RealExpression PID_setpoint(y=lambda_target)   annotation (Placement(transformation(extent={{-38,-2},{-18,18}})));
  Modelica.Blocks.Sources.BooleanExpression PID_control(y=lambda_control) annotation (Placement(transformation(extent={{-38,-20},{-18,0}})));
  ClaRa.Components.Utilities.Blocks.LimPID m_flow_PID(
    y_max=1,
    y_min=m_flow_rampup,
    Ni=N_i,
    y_inactive=0,
    use_activateInput=true,
    sign=1,
    Tau_i=tau_i,
    k=k_p,
    controllerType=Modelica.Blocks.Types.SimpleController.PI) annotation (Placement(transformation(extent={{2,-12},{26,12}})));

equation

  lambda_control = u1 > 0;

  // _____________________________________________
  //
  //           Characteristic Equations
  // _____________________________________________

  connect(PID_control.y, m_flow_PID.activateInput) annotation (Line(points={{-17,-10},{-10,-10},{-10,-9.6},{-0.4,-9.6}}, color={255,0,255}));
  connect(PID_setpoint.y, m_flow_PID.u_s) annotation (Line(points={{-17,8},{-6,8},{-6,0},{-0.4,0}}, color={0,0,127}));
  connect(u1, m_flow_PID.u_m) annotation (Line(points={{-102,-36},{14.12,-36},{14.12,-14.4}}, color={0,0,127}));
  connect(m_flow_PID.y, y) annotation (Line(points={{27.2,0},{108,0}}, color={0,0,127}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})),
                Documentation(info="<html>
<h4><span style=\"color: #008000\">1. Purpose of model</span></h4>
<p>Controller for Lambda in Fuel Cell Applications.</p>
<h4><span style=\"color: #008000\">2. Level of detail, physical effects considered, and physical insight</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">3. Limits of validity </span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">4. Interfaces</span></h4>
<p>Modelica RealInput: u1</p>
<p>Modelica RealOutput: y</p>
<h4><span style=\"color: #008000\">5. Nomenclature</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">6. Governing Equations</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">7. Remarks for Usage</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">8. Validation</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">9. References</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">10. Version History</span></h4>
<p>Model created by Simon Weilbach (simon.weilbach@tuhh.de) on 01.10.2014</p>
<p>Model revised by Pascal Dubucq (dubucq@tuhh.de) on 01.10.2015</p>
<p>Quality check (Code conventions) by Rebekka Denninger on 01.10.2016</p>
</html>"));
end LambdaController_PID;
