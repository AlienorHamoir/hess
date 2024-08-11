within H2Microgrid_TransiEnt.HESS;
model HESS_nonCompressed "HESS with low-pressure non-compressed storage system"

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;

  Modelica.Units.SI.Power P_FC_set "FC system electrical setpoint";
  Modelica.Units.SI.Power P_EL_set "Electrolyzer system electrical setpoint";
  Modelica.Blocks.Interfaces.RealOutput socTankH2 annotation (Placement(transformation(extent={{92,-18},{126,16}}), iconTransformation(extent={{92,-18},{126,16}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn P_set_HESS "Input for HESS power production and consumption setpoint" annotation (Placement(transformation(extent={{-132,-50},{-94,-12}}), iconTransformation(extent={{-132,-50},{-94,-12}})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_HESS "Actual HESS electrical power balance (consumed - produced)" annotation (Placement(transformation(extent={{92,60},{124,92}}), iconTransformation(extent={{92,60},{124,92}})));
  StorageSystem.H2StorageSystem_nonCompressed h2StorageSystem_nonCompressed annotation (Placement(transformation(extent={{20,-50},{56,6}})));
  ElectrolyzerBoPSystem.Electrolyzer.Systems.SystemElectrolyzerL2_nonCompressedStorage systemElectrolyzerL2_nonCompressedStorage(usePowerPort=false)
                                                                                                                                 annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-20,-60})));
  inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{64,-90},{84,-70}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_Txim_flow H2massSink(medium=medium, variable_m_flow=true) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={40,44})));
  Modelica.Blocks.Interfaces.RealInput T_environment "Prescribed boundary temperature from weather file" annotation (Placement(transformation(extent={{-134,2},{-94,42}})));
  Modelica.Blocks.Sources.RealExpression P_FC(y=P_FC_set) "Applied FC system electrical setpoint" annotation (Placement(transformation(extent={{-88,54},{-68,74}})));
  Modelica.Blocks.Sources.RealExpression P_EL(y=P_EL_set) "Applied electrolyzer system electrical setpoint" annotation (Placement(transformation(extent={{-86,-70},{-66,-50}})));
  Modelica.Blocks.Math.Add add annotation (Placement(transformation(extent={{50,66},{70,86}})));
equation

   if P_set_HESS >= 0 then
    P_EL_set = P_set_HESS;
    P_FC_set = 0;
  else
    P_EL_set = 0;
    P_FC_set = - P_set_HESS;
   end if;

  connect(systemElectrolyzerL2_nonCompressedStorage.gasPortOut, h2StorageSystem_nonCompressed.H2PortIn) annotation (Line(
      points={{-0.2,-60},{14,-60},{14,-22},{20,-22}},
      color={255,255,0},
      thickness=1.5));
  connect(h2StorageSystem_nonCompressed.socTank, socTankH2) annotation (Line(points={{57.08,-0.72},{92,-0.72},{92,-1},{109,-1}},
                                                                                                                               color={0,0,127}));
  connect(systemPEMFC.mflowH2_FC, H2massSink.m_flow) annotation (Line(points={{2,48.8},{20,48.8},{20,50},{28,50}},   color={0,0,127}));
  connect(H2massSink.gasPort, h2StorageSystem_nonCompressed.H2PortOut) annotation (Line(
      points={{50,44},{62,44},{62,-22},{56,-22}},
      color={255,255,0},
      thickness=1.5));
  connect(T_environment, systemPEMFC.T_environment) annotation (Line(points={{-114,22},{-68,22},{-68,50},{-50,50},{-50,77.6},{-41.6,77.6}},
                                                                                                                                      color={0,0,127}));
  connect(T_environment, systemElectrolyzerL2_nonCompressedStorage.T_environment) annotation (Line(points={{-114,22},{-68,22},{-68,-43.4},{-40.6,-43.4}}, color={0,0,127}));
  connect(P_FC.y, systemPEMFC.P_el_set) annotation (Line(points={{-67,64},{-66,60},{-41.6,60}},     color={0,0,127}));
  connect(P_EL.y, systemElectrolyzerL2_nonCompressedStorage.P_el_set) annotation (Line(points={{-65,-60},{-40.8,-60}}, color={0,0,127}));
  connect(systemPEMFC.P_FC_tot, add.u1) annotation (Line(
      points={{-12,80.8},{-12,88},{38,88},{38,82},{48,82}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(systemElectrolyzerL2_nonCompressedStorage.electrolyzerPowerOut, add.u2) annotation (Line(
      points={{-4.8,-39.2},{10,-39.2},{10,70},{48,70}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(add.y, P_HESS) annotation (Line(points={{71,76},{108,76}},                 color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={255,213,170},
          fillPattern=FillPattern.Solid), Text(
          extent={{-90,42},{90,-36}},
          textColor={0,0,0},
          textString="HESS
30 bar")}),                                                      Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<h4>1. Purpose of model</h4>
<p>Hydrogen Energy Storage System for microgrid applications, containing a fuel cell system, an electrolyzer system and low-pressure non-compressed storage system.</p>
<h4>2. Level of detail, physical effects considered, and physical insight</h4>
<p>HESS power setpoints maagement between fuel cell ad electrolyzer systems.</p>
<p>Link between fuel cell system (ideal gas) and storage system (real gas).</p>
<h4>3. Limits of validity </h4>
<h4>4. Interfaces</h4>
<h4>5. Nomenclature</h4>
<h4>(no remarks)</h4>
<h4>6. Governing Equations</h4>
<h4>7. Remarks for Usage</h4>
<h4>8. Validation</h4>
<p>Tested in &quot;H2Microgrid_TransiEnt.HESS.TestHESS&quot;</p>
<h4>9. References</h4>
<h4>10. Version History</h4>
<p>Model created by Ali&eacute;nor Hamoir in June 2024</p>
</html>"));
end HESS_nonCompressed;
