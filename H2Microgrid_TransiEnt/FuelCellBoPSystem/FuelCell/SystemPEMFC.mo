within H2Microgrid_TransiEnt.FuelCellBoPSystem.FuelCell;
model SystemPEMFC "Fuel cell system, with power controller"

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

  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{70,-98},{90,-78}})));
  TransiEnt.Components.Electrical.FuelCellSystems.FuelCell.PEM FC(
    T_nom=75 + 273,
    A_cell=0.0625,
    i_0=0.08,
    i_L=6e4,
    Ri=6e-5,
    cp=850,
    ka=0.3,
    i_Loss=0.2,
    no_Cells=20,
    usePowerPort=true,
    T_stack(start=80 + 273),
    v_n=0.733,
    is_Shutdown(start=true),
    redeclare TransiEnt.Components.Boundaries.Electrical.ApparentPower.PowerVoltage powerBoundary(Use_input_connector_v=false, v_boundary=PEM.v_n) "PV-Boundary for ApparentPowerPort",
    redeclare TransiEnt.Basics.Interfaces.Electrical.ActivePowerPort epp)   annotation (Placement(transformation(extent={{-20,-28},{22,14}})));
TransiEnt.Components.Electrical.FuelCellSystems.FuelCell.Controller.PowerController PowerController(k=1000) annotation (Placement(transformation(rotation=0, extent={{12,64},{-8,84}})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi AirSink(
    variable_p=false,
    variable_xi=false,
    p_const=1e5,
    T_const=200 + 273.15,
    medium=FC.Air)        annotation (Placement(transformation(
        extent={{-6,-7},{6,7}},
        rotation=180,
        origin={56,-41})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi SyngasSink(
    variable_p=false,
    variable_xi=false,
    p_const=1e5,
    T_const=200 + 273.15,
    medium=FC.Syngas) annotation (Placement(transformation(
        extent={{-6,-9},{6,9}},
        rotation=180,
        origin={56,45})));
TransiEnt.Components.Electrical.FuelCellSystems.FuelCell.Controller.LambdaController LambdaHController(
    T=0.01,
    m_flow_rampup=0.00025,
    k=0.002) annotation (Placement(transformation(extent={{20,-78},{0,-58}})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow SyngasSource(
    variable_T=false,
    m_flow_const=5.1e-2,
    variable_m_flow=true,
    variable_xi=false,
    xi_const={0.11,0.04,0.13,0.25,0.04,0.12},
    T_const=80 + 273.15,
    medium=FC.Syngas) annotation (Placement(transformation(extent={{-56,-3},{-40,14}})));
  TransiEnt.Basics.Interfaces.Gas.IdealGasTempPortOut drainh1 annotation (Placement(transformation(extent={{88,4},{108,24}})));
  Modelica.Blocks.Interfaces.RealOutput mflowH2_FC_set annotation (Placement(transformation(
        extent={{-11,-11},{11,11}},
        rotation=-90,
        origin={-19,-111}), iconTransformation(
        extent={{-11,-11},{11,11}},
        rotation=-90,
        origin={-19,-111})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn P_el_set "Input for power production setpoint" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={24,108})));
  replaceable TransiEnt.Basics.Interfaces.Electrical.ActivePowerPort epp if usePowerPort constrainedby TransiEnt.Basics.Interfaces.Electrical.PartialPowerPort "Choice of power port" annotation (choicesAllMatching=true, Dialog(tab="General", group="General"), Placement(transformation(extent={{-110,42},{-90,62}})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow AirSource(
    variable_m_flow=false,
    variable_xi=false,
    m_flow_const=2.55e-7,
    xi_const={0.01,0.7},
    T_const=25 + 273,
    medium=FC.Air)    annotation (Placement(transformation(
        extent={{6.5,-9},{-6.5,9}},
        rotation=180,
        origin={-53.5,-41})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_FC_actual annotation (Placement(transformation(extent={{96,-62},{116,-42}})));
equation

  connect(FC.lambda_H,LambdaHController. u1) annotation (Line(
      points={{22,-23.8},{30,-23.8},{30,-62},{20,-62}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(FC.feedh, SyngasSource.gas_a) annotation (Line(
      points={{-20,5.6},{-22,5.5},{-40,5.5}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.drainh, SyngasSink.gas_a) annotation (Line(
      points={{22,5.6},{44,5.6},{44,45},{50,45}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.draina, AirSink.gas_a) annotation (Line(
      points={{22,-19.6},{40,-19.6},{40,-41},{50,-41}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.feeda, AirSource.gas_a) annotation (Line(
      points={{-20,-19.6},{-44,-19.6},{-44,-41},{-47,-41}},
      color={255,170,85},
      thickness=0.5));
  connect(PowerController.V_cell, FC.V_stack) annotation (Line(points={{11,68.6},{66,68.6},{66,-7},{22,-7}},
                                                                                                    color={0,0,127}));
  connect(PowerController.y, FC.I_load) annotation (Line(points={{-9,74},{-62,74},{-62,-8.26},{-16.22,-8.26}},             color={0,0,127}));

  connect(FC.drainh, drainh1) annotation (Line(
      points={{22,5.6},{82,5.6},{82,14},{98,14}},
      color={255,170,85},
      thickness=1.5));
  connect(LambdaHController.y, mflowH2_FC_set) annotation (Line(points={{-0.8,-68},{-19,-68},{-19,-111}}, color={0,0,127}));
  connect(LambdaHController.y, SyngasSource.m_flow) annotation (Line(points={{-0.8,-68},{-74,-68},{-74,10},{-60,10},{-60,10.6},{-56,10.6}}, color={0,0,127}));
  connect(PowerController.deltaP, P_el_set) annotation (Line(points={{11,80},{24,80},{24,108}}, color={0,127,127}));
  connect(epp, FC.epp) annotation (Line(
      points={{-100,52},{1,52},{1,5.18}},
      color={0,135,135},
      thickness=0.5));
  connect(FC.P_el, P_FC_actual) annotation (Line(
      points={{-7.4,-28},{-8,-28},{-8,-52},{106,-52}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
    experiment(StopTime=40),
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
</html>"),
    Icon(graphics={                                                       Rectangle(
          extent={{-100,100},{102,-100}},
          lineColor={0,0,0},
          fillColor={162,29,33},
          fillPattern=FillPattern.Solid), Text(
          extent={{-74,42},{78,-38}},
          textColor={255,255,255},
          textStyle={TextStyle.Bold},
          textString="FC")}));
end SystemPEMFC;
