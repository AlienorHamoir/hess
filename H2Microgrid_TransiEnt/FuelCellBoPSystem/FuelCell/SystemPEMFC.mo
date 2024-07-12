within H2Microgrid_TransiEnt.FuelCellBoPSystem.FuelCell;
model SystemPEMFC "Fuel cell system, with auxiliaries and power controller"


  parameter TransiEnt.Basics.Media.Gases.Gas_VDIWA_SG7_var Syngas=TransiEnt.Basics.Media.Gases.Gas_VDIWA_SG7_var() "Medium model H2" annotation (Dialog(group="Fundamental Definitions"));
  parameter TransiEnt.Basics.Media.Gases.Gas_MoistAir Air=TransiEnt.Basics.Media.Gases.Gas_MoistAir() "Medium model of air" annotation (Dialog(group="Fundamental Definitions"));

  final parameter Modelica.Units.SI.SpecificEnergy NCV_H2 = 1.19951e8 "J/kg - Net calorific value of hydrogen at 25 C and 1 bar";
  final parameter Modelica.Units.SI.SpecificEnergy GCV_H2 = 1.41788e8 "J/kg - Gross calorific value of hydrogen at 25 C and 1 bar";
  final parameter Modelica.Units.SI.SpecificEnthalpy h0 = 356955 "J/kg - Specific enthalpy at 25 C and 1 bar";

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
        origin={-56,-110}), iconTransformation(
        extent={{-12,-12},{12,12}},
        rotation=-90,
        origin={-56,-110})));
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
    T_const=23.5 + 273.15,
    xi_const={0,0,0,0,1,0})
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
  PEMFC          FC(
    Syngas=TransiEnt.Basics.Media.Gases.Gas_VDIWA_SG7_var(),
    I_shutdown=10,
    T_nom(displayUnit="K") = 273.15 + 75,
    T_stack_max(displayUnit="K") = 273.15 + 80,
    T_cool_set(displayUnit="K") = 273.15 + 75,
    usePowerPort=false,
    useHeatPort=true,
    tau_i=0.1,
    k_p=1300)                     annotation (Placement(transformation(extent={{-14,-18},{18,12}})));
  Controller.LambdaController_PID          lambdaHController_PID(lambda_target=1.5, m_flow_rampup=1e-6)    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-76,-58})));
  Controller.LambdaController_PID          lambdaOController_PID(lambda_target=2.05, m_flow_rampup=1e-8)
                                                                 "Controller that outputs the required air mass flow rate to meet OER (oxygen excess ratio) target " annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-22,-42})));
  CoolingSystem.HeatPortCooling.CoolingModel coolingModel(k_p=1000, tau_i=0.5)
                                                          annotation (Placement(transformation(extent={{44,40},{64,60}})));
  AirSupplySystem.AirCompressorSystem AirCompressorSystem annotation (Placement(transformation(extent={{28,-84},{48,-64}})));
  Modelica.Blocks.Sources.RealExpression P_el_out(y=P_el_sys) annotation (Placement(transformation(extent={{68,-50},{88,-30}})));

  Controller.PowerController powerController(i_max=120, P_min=500) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-52,68})));
  Modelica.Blocks.Math.Add3 SumPower(
    k1=+1,
    k2=+1,
    k3=+1) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-10,78})));
  Modelica.Blocks.Interfaces.RealInput T_environment "Prescribed boundary temperature from weather file" annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={80,110}), iconTransformation(
        extent={{-10,30},{10,10}},
        rotation=-90,
        origin={68,108})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{66,-90},{86,-70}})));
  Modelica.Blocks.Sources.RealExpression P_el_out_aux(y=P_el_aux) annotation (Placement(transformation(extent={{68,-68},{88,-48}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_FC_aux annotation (Placement(transformation(extent={{92,-70},{116,-46}}), iconTransformation(extent={{92,-52},{116,-28}})));
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
  H_flow_H2_GCV = (mflowH2_FC*(GCV_H2 + (FC.h_haus - FC.h_hein)));
  H_flow_H2_NCV = (mflowH2_FC*(NCV_H2 + (FC.h_haus - FC.h_hein)));

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

  connect(FC.feedh,SyngasSource. gas_a) annotation (Line(
      points={{-14,6},{-30,6},{-30,19},{-36,19}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.feeda,AirSource. gas_a) annotation (Line(
      points={{-14,-12},{-32,-12},{-32,-19},{-37,-19}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.drainh,SyngasSink. gas_a) annotation (Line(
      points={{18,6},{38,6},{38,17},{50,17}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.draina,AirSink. gas_a) annotation (Line(
      points={{18,-12},{34,-12},{34,-20},{50,-20}},
      color={255,170,85},
      thickness=0.5));
  connect(FC.lambda_H,lambdaHController_PID. u1) annotation (Line(points={{10.32,-18},{10.32,-54.4},{-65.8,-54.4}},
                                                                                                              color={0,0,127}));
  connect(lambdaHController_PID.y,SyngasSource. m_flow) annotation (Line(points={{-86.8,-58},{-92,-58},{-92,24},{-62,24},{-62,23.8},{-52,23.8}}, color={0,0,127}));
  connect(FC.lambda_O,lambdaOController_PID. u1) annotation (Line(points={{-5.68,-18},{-5.68,-28},{-6,-28},{-6,-38},{-10,-38},{-10,-38.4},{-11.8,-38.4}},
                                                                                                                                                      color={0,0,127}));
  connect(FC.heat,coolingModel. heatPortCooling) annotation (Line(points={{18.16,-7.65},{66,-7.65},{66,36},{44,36},{44,42.2}},
                                                                                                                         color={191,0,0}));
  connect(FC.temperatureOut,coolingModel. T_op) annotation (Line(points={{-5.04,-7.5},{-22,-7.5},{-22,20},{36,20},{36,56},{44,56}},         color={0,0,127}));
  connect(lambdaOController_PID.y,AirSource. m_flow) annotation (Line(points={{-32.8,-42},{-56,-42},{-56,-24.4},{-50,-24.4}}, color={0,0,127}));
  connect(lambdaOController_PID.y, AirCompressorSystem.AirMassFlowRateSetpoint) annotation (Line(points={{-32.8,-42},{-36,-42},{-36,-67.6},{27.2,-67.6}}, color={0,0,127}));
  connect(lambdaHController_PID.y, mflowH2_FC) annotation (Line(points={{-86.8,-58},{-92,-58},{-92,-92},{-56,-92},{-56,-110}}, color={0,0,127}));
  connect(P_el_out.y, P_FC_tot) annotation (Line(points={{89,-40},{104,-40}}, color={0,0,127}));
  connect(FC.V_stack, powerController.V_stack) annotation (Line(
      points={{18,-3},{18,-4},{26,-4},{26,73.4},{-43,73.4}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(powerController.y, FC.I_load) annotation (Line(points={{-63,68},{-68,68},{-68,-3.9},{-11.12,-3.9}}, color={0,0,127}));
  connect(P_el_set, SumPower.u3) annotation (Line(points={{24,108},{24,86},{2,86}}, color={0,127,127}));
  connect(SumPower.y, powerController.P) annotation (Line(points={{-21,78},{-32,78},{-32,62},{-43,62}}, color={0,0,127}));
  connect(coolingModel.P_coolingPump, SumPower.u2) annotation (Line(points={{64,54},{68,54},{68,74},{28,74},{28,78},{2,78}},       color={0,0,127}));
  connect(AirCompressorSystem.P_airCompressor, SumPower.u1) annotation (Line(
      points={{48.6,-74.8},{56,-74.8},{56,-32},{44,-32},{44,22},{10,22},{10,70},{2,70}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(coolingModel.T_environment, T_environment) annotation (Line(points={{44,50},{34,50},{34,86},{80,86},{80,110}}, color={0,0,127}));
  connect(P_el_out_aux.y, P_FC_aux) annotation (Line(points={{89,-58},{104,-58}}, color={0,0,127}));
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
end SystemPEMFC;
