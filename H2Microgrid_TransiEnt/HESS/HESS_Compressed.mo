within H2Microgrid_TransiEnt.HESS;
model HESS_Compressed "HESS with high-pressure compressed storage system"

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;
  parameter Modelica.Units.SI.Pressure p_max=350e5   "Pressure in storage at t=0" annotation (Dialog(tab="General", group="Storage parameters"));
  parameter Real SOC_start=0.5  "Pressure in storage at t=0" annotation (Dialog(tab="General", group="Storage parameters"));
  parameter Modelica.Units.SI.Pressure p_start = SOC_start * p_max "Start pressure depending on inital HESS SOC";
  parameter Modelica.Units.SI.Power P_STB_FC = 16.6 "Constant power consumed by FC system due to auxiliaries in STANDBY";
  parameter Modelica.Units.SI.Power P_STB_EL = 285.7 "Constant power consumed by ELECTROLYZER system due to auxiliaries in STANDBY";

  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{28,82},{42,96}})));
  Modelica.Units.SI.Power P_cmd_FC(start=0) "FC system electrical setpoint";
  Modelica.Units.SI.Power P_cmd_EL(start=0) "Electrolyzer system electrical setpoint";
  Modelica.Units.SI.Power P_real_FC(start=0) "FC system actual electricity balance";
  Modelica.Units.SI.Power P_real_EL(start=0) "Electrolyzer system actual electricity balance";
  ElectrolyzerBoPSystem.Electrolyzer.Systems.SystemElectrolyzerL2_CompressedStorage systemEL(usePowerPort=false, electrolyzer(usePowerPort=false, integrateElPower=false)) annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-20,-48})));
  StorageSystem.H2StorageSystem_Compressed h2StorageSystem_Compressed(p_start=p_start, p_maxHigh=p_max)
                                                                      annotation (Placement(transformation(extent={{22,-40},{62,20}})));
  Modelica.Blocks.Interfaces.RealOutput LOH "H2 Tank Level of Hydrogen [-]"
                                                                          annotation (Placement(transformation(extent={{88,-16},{120,16}}), iconTransformation(extent={{88,-16},{120,16}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn P_set_FC "Input for FC power production setpoint" annotation (Placement(transformation(extent={{-120,50},{-82,90}}), iconTransformation(extent={{-120,50},{-82,90}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_Txim_flow H2massSink(medium=medium, variable_m_flow=true,
    T_const=313.15)                                                                                             annotation (Placement(transformation(
        extent={{-9,-8},{9,8}},
        rotation=0,
        origin={27,44})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_FC "Actual FC electrical power balance (consumed - produced)" annotation (Placement(transformation(extent={{90,48},{120,78}}), iconTransformation(extent={{90,48},{120,78}})));
  Modelica.Blocks.Interfaces.RealInput T_environment "Prescribed boundary temperature from weather file" annotation (Placement(transformation(extent={{-114,-20},{-74,20}}), iconTransformation(extent={{-114,-20},{-74,20}})));
  inner TransiEnt.ModelStatistics modelStatistics annotation (Placement(transformation(extent={{48,82},{62,96}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn P_set_EL "Input for Electrolyzer power consumption setpoint" annotation (Placement(transformation(extent={{-120,-90},{-82,-50}}), iconTransformation(extent={{-120,-90},{-82,-50}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_EL "Actual Electrolyzer electrical power balance (consumed - produced)" annotation (Placement(transformation(extent={{92,-82},{122,-52}}), iconTransformation(extent={{92,-82},{122,-52}})));
  Modelica.Blocks.Sources.RealExpression PowerSetpointFC(y=P_cmd_FC) "Applied FC system electrical setpoint" annotation (Placement(transformation(extent={{-70,50},{-50,70}})));
  Modelica.Blocks.Sources.RealExpression PowerSetpointEL(y=P_cmd_EL) "Applied electrolyzer system electrical setpoint" annotation (Placement(transformation(extent={{-76,-58},{-56,-38}})));
  Modelica.Blocks.Sources.RealExpression PowerOutEL(y=P_real_EL) "Actual electricity consumed by electrolyzer system" annotation (Placement(transformation(extent={{42,-76},{62,-56}})));
  Modelica.Blocks.Sources.RealExpression PowerOutFC(y=P_real_FC) "Actual electricity produced/consumed by FC system"
                                                                                                            annotation (Placement(transformation(extent={{56,52},{76,72}})));
  Modelica.Blocks.Interfaces.RealInput state_FC "State of the fuel cell (0 = OFF, 1 = STANDBY, 2 = ON, 3 = ERROR)" annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-40,-100})));
  Modelica.Blocks.Interfaces.RealInput state_EL "State of the electrolyzer (0 = OFF, 1 = STANDBY, 2 = ON, 3 = ERROR)" annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={40,-100})));
  FuelCellBoPSystem.FuelCell.SystemPEMFCexp systemFC annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-20,60})));
equation
  //// States and power setpoints management: delta_i = 0 -> OFF, delta_i = 1 -> STB, delta_i = 2 -> ON
  // FUEL CELL
  if state_FC == 0 then
    P_cmd_FC = 0;
    P_real_FC = 0;
  elseif state_FC == 1 then
    P_cmd_FC = 0;
    P_real_FC = P_STB_FC;
  else
    P_cmd_FC = P_set_FC;
    P_real_FC =systemFC.P_FC_sys;
  end if;

  // ELECTROLYZER
  if state_EL == 0 then
    P_cmd_EL = 0;
    P_real_EL = 0;
  elseif state_EL == 1 then
    P_cmd_EL = 0;
    P_real_EL = P_STB_EL;
  else
    P_cmd_EL = P_set_EL;
    P_real_EL = systemEL.P_electrolyzer_tot;
  end if;

  connect(systemEL.gasPortOut, h2StorageSystem_Compressed.H2PortIn) annotation (Line(
      points={{-0.2,-48},{41.8,-48},{41.8,-40.9}},
      color={255,255,0},
      thickness=1.5));
  connect(h2StorageSystem_Compressed.LOH, LOH) annotation (Line(points={{62,13.4},{70,13.4},{70,14},{78,14},{78,0},{104,0}},   color={0,0,127}));
  connect(h2StorageSystem_Compressed.H2PortOut, H2massSink.gasPort) annotation (Line(
      points={{42,20.6},{42,44},{36,44}},
      color={255,255,0},
      thickness=1.5));
  connect(h2StorageSystem_Compressed.P_comp, systemEL.CompressorPower) annotation (Line(
      points={{62,-20.8},{68,-20.8},{68,-76},{-7.8,-76},{-7.8,-68.6}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(systemEL.T_environment, T_environment) annotation (Line(points={{-40.6,-31.4},{-74,-31.4},{-74,0},{-94,0}},  color={0,0,127}));
  connect(PowerSetpointEL.y, systemEL.P_el_set) annotation (Line(points={{-55,-48},{-40.8,-48}},                     color={0,0,127}));
  connect(PowerOutEL.y, P_EL) annotation (Line(points={{63,-66},{88,-66},{88,-67},{107,-67}}, color={0,0,127}));
  connect(PowerOutFC.y, P_FC) annotation (Line(points={{77,62},{90,62},{90,63},{105,63}},
                                                                                  color={0,0,127}));
  connect(PowerSetpointFC.y,systemFC.P_FC_set)  annotation (Line(points={{-49,60},{-41.6,60}}, color={0,0,127}));
  connect(T_environment, systemFC.T_env) annotation (Line(points={{-94,0},{-74,0},{-74,84},{-41.6,84},{-41.6,77.6}}, color={0,0,127}));
  connect(systemFC.mflowH2, H2massSink.m_flow) annotation (Line(points={{0,48.8},{16.2,48.8}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={255,213,170},
          fillPattern=FillPattern.Solid), Text(
          extent={{-90,42},{90,-36}},
          textColor={0,0,0},
          textString="HESS
HP")}),                                                          Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<h4>1. Purpose of model</h4>
<p>Hydrogen Energy Storage System for microgrid applications, containing a fuel cell system, an electrolyzer system and high-pressure compressed storage system.</p>
<h4>2. Level of detail, physical effects considered, and physical insight</h4>
<p>HESS power setpoints maagement betwee fuel cell ad electrolyzer systems.</p>
<p>Link between fuel cell system (ideal gas) and storage system (real gas).</p>
<h4>3. Limits of validity </h4>
<h4>4. Interfaces</h4>
<h4>5. Nomenclature</h4>
<p>(no remarks)</p>
<h4>6. Governing Equations</h4>
<h4>7. Remarks for Usage</h4>
<h4>8. Validation</h4>
<p>Tested in &quot;H2Microgrid_TransiEnt.HESS.TestHESS&quot;</p>
<h4>9. References</h4>
<h4>10. Version History</h4>
<p>Model created by Ali&eacute;nor Hamoir in June 2024</p>
</html>"));
end HESS_Compressed;
