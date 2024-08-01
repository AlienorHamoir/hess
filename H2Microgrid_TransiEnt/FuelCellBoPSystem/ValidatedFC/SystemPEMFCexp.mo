within H2Microgrid_TransiEnt.FuelCellBoPSystem.ValidatedFC;
model SystemPEMFCexp "Fuel cell system, with auxiliaries and power controller - experimentally validated fuel cell model"

  parameter TransiEnt.Basics.Media.Gases.Gas_VDIWA_SG7_var Syngas=TransiEnt.Basics.Media.Gases.Gas_VDIWA_SG7_var() "Medium model H2" annotation (Dialog(group="Fundamental Definitions"));
  parameter TransiEnt.Basics.Media.Gases.Gas_MoistAir Air=TransiEnt.Basics.Media.Gases.Gas_MoistAir() "Medium model of air" annotation (Dialog(group="Fundamental Definitions"));

  final parameter Modelica.Units.SI.SpecificEnergy NCV_H2 = 1.19951e8 "J/kg - Net calorific value of hydrogen at 25 C and 1 bar";
  final parameter Modelica.Units.SI.SpecificEnergy GCV_H2 = 1.41788e8 "J/kg - Gross calorific value of hydrogen at 25 C and 1 bar";
  final parameter Modelica.Units.SI.SpecificEnthalpy h0 = 356955 "J/kg - Reference specific enthalpy at 25 C and 1 bar";

  Modelica.Units.SI.Power P_el_tot "Electric power consumed by the fuel cell and air compressor and cooling";
  Modelica.Units.SI.Power P_el_sys "Electric power actually consumed by the system FC and air compressor and cooling ";
  Modelica.Units.SI.Power P_el_aux "Electric power consumed by the fuel cell auxiliaries system";
  Modelica.Units.SI.EnthalpyFlowRate H_flow_H2_GCV "H2 enthalpy flow rate out of electrolyzer, gross calorific value";
  Modelica.Units.SI.EnthalpyFlowRate H_flow_H2_NCV "H2 enthalpy flow rate out of electrolyzer, net calorific value";
  Modelica.Units.SI.Efficiency eta_NCV_FC(min=0, max=1) "Efficiency of the electrolyzer based on NCV" annotation (Dialog(group="Initialization", showStartAttribute=true));
  Modelica.Units.SI.Efficiency eta_GCV_FC(min=0, max=1) "Efficiency of the electrolyzer based on GCV" annotation (Dialog(group="Initialization", showStartAttribute=true));
  Modelica.Units.SI.Efficiency eta_NCV_sys(min=0, max=1) "Efficiency of the electrolyzer + BoP system based on NCV" annotation (Dialog(group="Initialization", showStartAttribute=true));
  Modelica.Units.SI.Efficiency eta_GCV_sys(min=0, max=1) "Efficiency of the electrolyzer + BoP system based on GCV" annotation (Dialog(group="Initialization", showStartAttribute=true));

  Modelica.Blocks.Interfaces.RealOutput mflowH2_FC annotation (Placement(transformation(
        extent={{-12,-12},{12,12}},
        rotation=-90,
        origin={-56,-102}), iconTransformation(
        extent={{-12,-12},{12,12}},
        rotation=-90,
        origin={-56,-102})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn P_el_set "Input for FC system power production setpoint" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={24,108}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={0,108})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_FC_tot annotation (Placement(transformation(extent={{92,-52},{116,-28}}), iconTransformation(extent={{92,-52},{116,-28}})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow AirSource(
    variable_m_flow=true,
    variable_xi=false,
    m_flow_const=0.001,
    T_const=23.5 + 273,
    medium=FC.Air)    annotation (Placement(transformation(
        extent={{6.5,-9},{-6.5,9}},
        rotation=180,
        origin={-43.5,-19})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi AirSink(
    variable_p=false,
    variable_T=false,
    variable_xi=false,
    p_const=1e5,
    T_const=23.50 + 273.15,
    medium=FC.Air)        annotation (Placement(transformation(
        extent={{-7,-8},{7,8}},
        rotation=180,
        origin={57,-20})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi SyngasSink(
    variable_p=false,
    variable_T=false,
    variable_xi=false,
    p_const=1e5,
    medium=FC.Syngas,
    T_const=23.5 + 273.15)
                      annotation (Placement(transformation(
        extent={{-6,-9},{6,9}},
        rotation=180,
        origin={56,17})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow SyngasSource(
    variable_T=false,
    m_flow_const=5.1e-2,
    variable_m_flow=true,
    variable_xi=false,
    T_const=40 + 273.15,
    medium=FC.Syngas,
    xi_const={0,0,0,0,1,0})
                      annotation (Placement(transformation(extent={{-52,11},{-36,27}})));
  FuelCell.Controller.LambdaController_PID lambdaOController_PID(lambda_target=2,    m_flow_rampup=1e-8) "Controller that outputs the required air mass flow rate to meet OER (oxygen excess ratio) target " annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-22,-42})));
  CoolingSystem.HeatPortCooling.CoolingModel coolingModel(k_p=1000, tau_i=0.01)
                                                          annotation (Placement(transformation(extent={{52,46},{72,66}})));
  AirSupplySystem.AirCompressorSystem AirCompressorSystem annotation (Placement(transformation(extent={{16,-74},{36,-54}})));
  Modelica.Blocks.Sources.RealExpression P_el_out(y=P_el_sys) annotation (Placement(transformation(extent={{68,-50},{88,-30}})));

  FuelCell.Controller.PowerController powerController(i_max=200, preventDivisionByZero(uMax=45))
                                                                 annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-38,68})));
  Modelica.Blocks.Interfaces.RealInput T_environment "Prescribed boundary temperature from weather file" annotation (Placement(transformation(
        extent={{-15,-15},{15,15}},
        rotation=-90,
        origin={85,107}), iconTransformation(
        extent={{-10,30},{10,10}},
        rotation=-90,
        origin={68,108})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{62,-98},{76,-84}})));
  Modelica.Blocks.Sources.RealExpression P_el_out_aux(y=P_el_aux) annotation (Placement(transformation(extent={{68,-88},{88,-68}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_FC_aux annotation (Placement(transformation(extent={{92,-92},{116,-68}}), iconTransformation(extent={{92,-92},{116,-68}})));
  PEMFC_KhanIqbal FC(useHeatPort=true) annotation (Placement(transformation(extent={{-12,-14},{16,12}})));
  FuelCell.Controller.LambdaController_PID lambdaHController_PID(lambda_target=1.5, m_flow_rampup=1e-8) "Controller that outputs the required air mass flow rate to meet OER (oxygen excess ratio) target " annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-76,-58})));
  Modelica.Blocks.Math.Add3 SumPower(
    k1=+1,
    k2=+1,
    k3=+1) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-2,76})));
equation

  P_el_tot = FC.P_el + AirCompressorSystem.P_airCompressor + coolingModel.P_coolingPump;
  P_el_aux = AirCompressorSystem.P_airCompressor + coolingModel.P_coolingPump;

  //System Power balance
  if abs(FC.P_el)>0 then
    P_el_sys = P_el_tot;
  else
    P_el_sys = 0;
  end if;

  //Efficiency
  H_flow_H2_GCV = (mflowH2_FC*(GCV_H2 + (h0 + FC.h_hein)));
  H_flow_H2_NCV = (mflowH2_FC*(NCV_H2 + (h0 + FC.h_hein)));

  if abs(FC.P_el)>0 then
    eta_GCV_FC = - FC.P_el / H_flow_H2_GCV;
    eta_NCV_FC = - FC.P_el / H_flow_H2_NCV;
    eta_GCV_sys = - P_el_tot / H_flow_H2_GCV;
    eta_NCV_sys = - P_el_tot / H_flow_H2_NCV;
  else
    eta_GCV_FC = 0;
    eta_NCV_FC = eta_GCV_FC;
    eta_GCV_sys = 0;
    eta_NCV_sys = eta_GCV_sys;
  end if;

  //// Discontinuous peak of efficiency higher than 1: can be managed by post processing in MPC
//   if eta_GCV_sys > 1
//     eta_NCV_FC = 1;
//     eta_NCV_sys = 1;
//     eta_GCV_FC = 1;
//     eta_GCV_sys = 1;
//   else
//     eta_NCV_FC = eta_NCV_FC;
//     eta_NCV_sys = eta_NCV_sys;
//     eta_GCV_FC = eta_GCV_FC;
//     eta_GCV_sys = eta_GCV_sys;
//   end if;

  connect(FC.feedh, SyngasSource.gas_a) annotation (Line(
      points={{-12,6.8},{-28,6.8},{-28,19},{-36,19}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.feeda, AirSource.gas_a) annotation (Line(
      points={{-12,-8.8},{-28,-8.8},{-28,-19},{-37,-19}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.drainh, SyngasSink.gas_a) annotation (Line(
      points={{16,6.8},{16,17},{50,17}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.draina, AirSink.gas_a) annotation (Line(
      points={{16,-8.8},{50,-8.8},{50,-20},{50,-20}},
      color={255,170,85},
      thickness=0.5));
  connect(lambdaOController_PID.y,AirSource. m_flow) annotation (Line(points={{-31.4,-42},{-56,-42},{-56,-24.4},{-50,-24.4}}, color={0,0,127}));
  connect(lambdaOController_PID.y, AirCompressorSystem.AirMassFlowRateSetpoint) annotation (Line(points={{-31.4,-42},{-36,-42},{-36,-58},{16,-58}},       color={0,0,127}));
  connect(P_el_out.y, P_FC_tot) annotation (Line(points={{89,-40},{104,-40}}, color={0,0,127}));
  connect(coolingModel.T_environment, T_environment) annotation (Line(points={{52,56},{34,56},{34,86},{85,86},{85,107}}, color={0,0,127}));
  connect(P_el_out_aux.y, P_FC_aux) annotation (Line(points={{89,-78},{90,-80},{104,-80}},          color={0,0,127}));
  connect(FC.lambda_O, lambdaOController_PID.u1) annotation (Line(points={{-4.72,-14},{-10,-14},{-10,-38},{-12.6,-38}}, color={0,0,127}));
  connect(powerController.y, FC.I_load) annotation (Line(points={{-49,68},{-58,68},{-58,-1.78},{-9.48,-1.78}},                   color={0,0,127}));
  connect(FC.V_stack, powerController.V_stack) annotation (Line(
      points={{16,-1},{16,-2},{28,-2},{28,66},{-22,66},{-22,73.4},{-29,73.4}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(FC.heat, coolingModel.heatPortCooling) annotation (Line(points={{16.14,-5.03},{46,-5.03},{46,48.2},{52,48.2}},                 color={191,0,0}));
  connect(FC.temperatureOut, coolingModel.T_op) annotation (Line(points={{-4.16,-4.9},{-18,-4.9},{-18,18},{26,18},{26,62},{52,62}}, color={0,0,127}));
  connect(FC.lambda_H, lambdaHController_PID.u1) annotation (Line(points={{9.28,-14},{9.28,-60},{-60,-60},{-60,-54},{-66.6,-54}}, color={0,0,127}));
  connect(lambdaHController_PID.y, SyngasSource.m_flow) annotation (Line(points={{-85.4,-58},{-90,-58},{-90,23.8},{-52,23.8}}, color={0,0,127}));
  connect(lambdaHController_PID.y, mflowH2_FC) annotation (Line(points={{-85.4,-58},{-90,-58},{-90,-84},{-56,-84},{-56,-102}}, color={0,0,127}));
  connect(SumPower.y, powerController.P) annotation (Line(points={{-13,76},{-20,76},{-20,62},{-29,62}}, color={0,0,127}));
  connect(SumPower.u1, AirCompressorSystem.P_airCompressor) annotation (Line(points={{10,68},{10,20},{-20,20},{-20,-24},{42,-24},{42,-66},{36,-66}}, color={0,0,127}));
  connect(SumPower.u2, coolingModel.P_coolingPump) annotation (Line(points={{10,76},{82,76},{82,60},{72,60}}, color={0,0,127}));
  connect(SumPower.u3, P_el_set) annotation (Line(points={{10,84},{10,94},{24,94},{24,108}}, color={0,0,127}));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
    experiment(
      StopTime=2500,
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
<p>PEMFC system model, with OER and HER controllers, power controller, air compressor system and cooling system models</p>
<h4><span style=\"color: #008000\">2. Level of detail, physical effects considered, and physical insight</span></h4>
<p> Efficiency equations</p>
<p>Power setpoint managemet</p>
<h4><span style=\"color: #008000\">3. Limits of validity </span></h4>
<p>(Purely technical component without physical modeling.)</p>
<h4><span style=\"color: #008000\">4. Interfaces</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">5. Nomenclature</span></h4>
<p>(no elements)</p>
<h4><span style=\"color: #008000\">6. Governing Equations</span></h4>
<p>NCV and GCV efficiencies equations </p>
<h4><span style=\"color: #008000\">7. Remarks for Usage</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">8. Validation</span></h4>
<p>(no validation or testing necessary)</p>
<h4><span style=\"color: #008000\">9. References</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">10. Version History</span></h4>
<p>Created by Ali&eacute;nor Hamoir in June 2024</p>
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
end SystemPEMFCexp;
