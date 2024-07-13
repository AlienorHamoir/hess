within H2Microgrid_TransiEnt.HESS;
model HESS_Compressed "HESS with high-pressure compressed storage system"

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;
  parameter Modelica.Units.SI.Pressure p_max=350e5   "Pressure in storage at t=0" annotation (Dialog(tab="General", group="Storage parameters"));
  parameter Real SOC_start=0.5  "Pressure in storage at t=0" annotation (Dialog(tab="General", group="Storage parameters"));
  parameter Modelica.Units.SI.Pressure p_start = SOC_start * p_max "Start pressure depending on inital HESS SOC";





inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{-92,-98},{-78,-84}})));
  Modelica.Units.SI.Power P_FC_set(start=0) "FC system electrical setpoint";
  Modelica.Units.SI.Power P_EL_set(start=0) "Electrolyzer system electrical setpoint";

  ElectrolyzerBoPSystem.Electrolyzer.Systems.SystemElectrolyzerL2_CompressedStorage systemElectrolyzer(usePowerPort=false, electrolyzer(usePowerPort=false, integrateElPower=false)) annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-20,-60})));
  StorageSystem.H2StorageSystem_Compressed h2StorageSystem_Compressed(p_start=p_start, p_maxHigh=p_max)
                                                                      annotation (Placement(transformation(extent={{22,-40},{62,20}})));
  Modelica.Blocks.Interfaces.RealOutput socTankH2 annotation (Placement(transformation(extent={{94,-42},{114,-22}})));
  FuelCellBoPSystem.FuelCell.SystemPEMFC systemPEMFC annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-20,60})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn P_set_HESS "Input for HESS power production or consumption setpoint" annotation (Placement(transformation(extent={{-124,-50},{-86,-10}}), iconTransformation(extent={{-124,-50},{-86,-10}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_Txim_flow H2massSink(medium=medium, variable_m_flow=true) annotation (Placement(transformation(
        extent={{-9,-8},{9,8}},
        rotation=0,
        origin={27,44})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_HESS "Actual HESS electrical power balance (consumed - produced)"
                                                                                                           annotation (Placement(transformation(extent={{96,52},{116,72}})));
  Modelica.Blocks.Math.Add add annotation (Placement(transformation(extent={{62,52},{82,72}})));
  Modelica.Blocks.Sources.RealExpression P_FC(y=P_FC_set) "Applied FC system electrical setpoint" annotation (Placement(transformation(extent={{-74,50},{-54,70}})));
  Modelica.Blocks.Sources.RealExpression P_EL(y=P_EL_set) "Applied electrolyzer system electrical setpoint" annotation (Placement(transformation(extent={{-72,-70},{-52,-50}})));
  Modelica.Blocks.Interfaces.RealInput T_environment "Prescribed boundary temperature from weather file" annotation (Placement(transformation(extent={{-124,10},{-84,50}})));
  inner TransiEnt.ModelStatistics modelStatistics annotation (Placement(transformation(extent={{-72,-98},{-58,-84}})));
equation

  if P_set_HESS >= 0 then
    P_EL_set = P_set_HESS;
    P_FC_set = 0;
  else
    P_EL_set = 0;
    P_FC_set = - P_set_HESS;
  end if;

  connect(systemElectrolyzer.gasPortOut, h2StorageSystem_Compressed.H2PortIn) annotation (Line(
      points={{-0.2,-60},{42.2,-60},{42.2,-40.3}},
      color={255,255,0},
      thickness=1.5));
  connect(h2StorageSystem_Compressed.socTank, socTankH2) annotation (Line(points={{63.2,13.4},{90,13.4},{90,-32},{104,-32}},
                                                                                                                         color={0,0,127}));
  connect(h2StorageSystem_Compressed.H2PortOut, H2massSink.gasPort) annotation (Line(
      points={{42,19.4},{42,44},{36,44}},
      color={255,255,0},
      thickness=1.5));
  connect(systemPEMFC.mflowH2_FC, H2massSink.m_flow) annotation (Line(points={{0.4,48.8},{16.2,48.8}},      color={0,0,127}));
  connect(h2StorageSystem_Compressed.P_comp, systemElectrolyzer.CompressorPower) annotation (Line(
      points={{62.8,-20.8},{70,-20.8},{70,-88},{-7.8,-88},{-7.8,-80.6}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(systemPEMFC.P_FC_tot, add.u1) annotation (Line(
      points={{-12,80.8},{-12,88},{52,88},{52,68},{60,68}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(systemElectrolyzer.P_electrolyzer_tot, add.u2) annotation (Line(
      points={{-12.4,-38.8},{-12.4,28},{52,28},{52,56},{60,56}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(add.y, P_HESS) annotation (Line(points={{83,62},{106,62}}, color={0,0,127}));
  connect(P_FC.y, systemPEMFC.P_el_set) annotation (Line(points={{-53,60},{-41.6,60}},              color={0,0,127}));
  connect(P_EL.y, systemElectrolyzer.P_el_set) annotation (Line(points={{-51,-60},{-40.8,-60}},                     color={0,0,127}));
  connect(systemElectrolyzer.T_environment, T_environment) annotation (Line(points={{-40.6,-43.4},{-80,-43.4},{-80,30},{-104,30}}, color={0,0,127}));
  connect(T_environment, systemPEMFC.T_environment) annotation (Line(points={{-104,30},{-80,30},{-80,77.6},{-41.6,77.6}},
                                                                                                                    color={0,0,127}));
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
