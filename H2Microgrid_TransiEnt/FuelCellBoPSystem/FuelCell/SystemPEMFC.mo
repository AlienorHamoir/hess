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

  parameter H2Microgrid_TransiEnt.FuelCellBoPSystem.FuelCell.Physics.Gas_VDIWA_H2_var Syngas=H2Microgrid_TransiEnt.FuelCellBoPSystem.FuelCell.Physics.Gas_VDIWA_H2_var() "Medium model H2" annotation (Dialog(group="Fundamental Definitions"));


  parameter TransiEnt.Basics.Media.Gases.Gas_MoistAir Air=TransiEnt.Basics.Media.Gases.Gas_MoistAir() "Medium model of air" annotation (Dialog(group="Fundamental Definitions"));

  Modelica.Blocks.Interfaces.RealOutput mflowH2_FC_set annotation (Placement(transformation(
        extent={{-12,-12},{12,12}},
        rotation=-90,
        origin={-36,-110}), iconTransformation(
        extent={{-12,-12},{12,12}},
        rotation=-90,
        origin={-36,-110})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn P_el_set "Input for power production setpoint" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={24,108})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_FC_actual annotation (Placement(transformation(extent={{100,-46},{120,-26}})));
TransiEnt.Components.Electrical.FuelCellSystems.FuelCell.Controller.PowerController PowerController(k=1, i_min=20)
                                                                                                            annotation (Placement(transformation(rotation=0, extent={{-38,62},{-58,82}})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow AirSource(
    variable_m_flow=true,
    variable_xi=false,
    m_flow_const=0.001,
    T_const=25 + 273,
    medium=FC.Air)    annotation (Placement(transformation(
        extent={{6.5,-9},{-6.5,9}},
        rotation=180,
        origin={-43.5,-11})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi AirSink(
    variable_p=false,
    variable_T=false,
    variable_xi=false,
    p_const=1e5,
    T_const=25 + 273.15,
    medium=FC.Air)        annotation (Placement(transformation(
        extent={{-7,-8},{7,8}},
        rotation=180,
        origin={57,-12})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi SyngasSink(
    variable_p=false,
    variable_T=false,
    variable_xi=false,
    p_const=1e5,
    medium=FC.Syngas,
    T_const=25 + 273.15,
    xi_const={0,0,0,0.005,0.99,0})
                      annotation (Placement(transformation(
        extent={{-6,-9},{6,9}},
        rotation=180,
        origin={56,25})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow SyngasSource(
    variable_T=false,
    m_flow_const=5.1e-2,
    variable_m_flow=true,
    variable_xi=false,
    T_const=40 + 273.15,
    medium=FC.Syngas,
    xi_const={0,0,0,0.005,0.99,0})
                      annotation (Placement(transformation(extent={{-52,19},{-36,35}})));
  PEMFC          FC(
    Syngas=TransiEnt.Basics.Media.Gases.Gas_VDIWA_SG7_var(),
    m=1,
    cp=35000,
    no_Cells=35,
    T_nom=308.15,
    T_stack_max=313.15,
    T_cool_set=303.15,
    T_stack(start=298.15),
    usePowerPort=false,
    useHeatPort=true)             annotation (Placement(transformation(extent={{-14,-10},{18,20}})));
  Controller.LambdaController_PID lambdaHController_PID(lambda_target=1.5, m_flow_rampup=1e-6) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-76,-50})));
  Controller.LambdaController_PID lambdaOController_PID annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-22,-34})));
  CoolingSystem.HeatPortCooling.CoolingModel coolingModel annotation (Placement(transformation(extent={{66,64},{86,84}})));
  AirSupplySystem.AirCompressor airCompressor annotation (Placement(transformation(
        extent={{-10,-9},{10,9}},
        rotation=-90,
        origin={-53,-68})));
equation

  connect(FC.feedh,SyngasSource. gas_a) annotation (Line(
      points={{-14,14},{-30,14},{-30,27},{-36,27}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.feeda,AirSource. gas_a) annotation (Line(
      points={{-14,-4},{-32,-4},{-32,-11},{-37,-11}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.drainh,SyngasSink. gas_a) annotation (Line(
      points={{18,14},{38,14},{38,25},{50,25}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.draina,AirSink. gas_a) annotation (Line(
      points={{18,-4},{34,-4},{34,-12},{50,-12}},
      color={255,170,85},
      thickness=0.5));
  connect(PowerController.V_cell,FC. V_stack) annotation (Line(points={{-39,66.6},{28,66.6},{28,5},{18,5}},
                                                                                                    color={0,0,127}));
  connect(PowerController.y,FC. I_load) annotation (Line(points={{-59,72},{-66,72},{-66,4},{-38,4},{-38,4.1},{-11.12,4.1}},color={0,0,127}));
  connect(FC.lambda_H,lambdaHController_PID. u1) annotation (Line(points={{10.32,-10},{10.32,-46.4},{-65.8,-46.4}},
                                                                                                              color={0,0,127}));
  connect(lambdaHController_PID.y,SyngasSource. m_flow) annotation (Line(points={{-86.8,-50},{-92,-50},{-92,32},{-62,32},{-62,31.8},{-52,31.8}}, color={0,0,127}));
  connect(FC.lambda_O,lambdaOController_PID. u1) annotation (Line(points={{-5.68,-10},{-5.68,-20},{-6,-20},{-6,-30},{-10,-30},{-10,-30.4},{-11.8,-30.4}},
                                                                                                                                                      color={0,0,127}));
  connect(FC.heat,coolingModel. heatPortCooling) annotation (Line(points={{18.16,0.35},{66,0.35},{66,64.8}},             color={191,0,0}));
  connect(FC.temperatureOut,coolingModel. T_op) annotation (Line(points={{-5.04,0.5},{6,0.5},{6,52},{60,52},{60,80.6},{64.8,80.6}},         color={0,0,127}));
  connect(lambdaOController_PID.y, airCompressor.AirMassFlowRateSetpoint) annotation (Line(points={{-32.8,-34},{-36,-34},{-36,-50},{-45.8,-50},{-45.8,-57}}, color={0,0,127}));
  connect(lambdaOController_PID.y, AirSource.m_flow) annotation (Line(points={{-32.8,-34},{-60,-34},{-60,-16.4},{-50,-16.4}}, color={0,0,127}));
  connect(P_el_set, PowerController.deltaP) annotation (Line(points={{24,108},{24,78},{-39,78}}, color={0,127,127}));
  connect(FC.P_el, P_FC_actual) annotation (Line(
      points={{18.96,-9.4},{18.96,-36},{110,-36}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(lambdaOController_PID.y, mflowH2_FC_set) annotation (Line(points={{-32.8,-34},{-36,-34},{-36,-110}}, color={0,0,127}));
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
