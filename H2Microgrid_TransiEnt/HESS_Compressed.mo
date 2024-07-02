within H2Microgrid_TransiEnt;
model HESS_Compressed

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.gasModel3;

inner TransiEnt.SimCenter simCenter annotation (Placement(transformation(extent={{-88,-98},{-68,-78}})));
  Modelica.Units.SI.Power P_FC_set "FC system electrical setpoint";
  Modelica.Units.SI.Power P_EL_set "Electrolyzer system electrical setpoint";

  ElectrolyzerBoPSystem.Electrolyzer.Systems.SystemElectrolyzerL2_CompressedStorage systemElectrolyzer(usePowerPort=false, electrolyzer(usePowerPort=false, integrateElPower=false)) annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-20,-60})));
  StorageSystem.H2StorageSystem_Compressed h2StorageSystem_Compressed annotation (Placement(transformation(extent={{22,-40},{62,20}})));
  Modelica.Blocks.Interfaces.RealOutput socTankH2 annotation (Placement(transformation(extent={{94,-42},{114,-22}})));
  FuelCellBoPSystem.FuelCell.SystemPEMFC systemPEMFC annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-20,60})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerIn P_set_HESS "Input for HESS power production or consumption setpoint" annotation (Placement(transformation(extent={{-128,-50},{-90,-10}}), iconTransformation(extent={{-128,-50},{-90,-10}})));
  TransiEnt.Components.Boundaries.Gas.BoundaryRealGas_Txim_flow H2massSink(medium=medium, variable_m_flow=true) annotation (Placement(transformation(
        extent={{-9,-8},{9,8}},
        rotation=0,
        origin={27,44})));
  TransiEnt.Basics.Interfaces.Electrical.ElectricPowerOut P_HESS "Actual HESS electrical power production" annotation (Placement(transformation(extent={{96,52},{116,72}})));
  Modelica.Blocks.Math.Add add annotation (Placement(transformation(extent={{62,52},{82,72}})));
  Modelica.Blocks.Sources.RealExpression P_FC(y=P_FC_set) "Applied FC system electrical setpoint" annotation (Placement(transformation(extent={{-76,54},{-56,74}})));
  Modelica.Blocks.Sources.RealExpression P_EL(y=P_EL_set) "Applied electrolyzer system electrical setpoint" annotation (Placement(transformation(extent={{-74,-70},{-54,-50}})));
  Modelica.Blocks.Interfaces.RealInput T_environment "Prescribed boundary temperature from weather file" annotation (Placement(transformation(extent={{-130,10},{-90,50}})));
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
  connect(systemPEMFC.mflowH2_FC, H2massSink.m_flow) annotation (Line(points={{2,48.8},{16.2,48.8}},        color={0,0,127}));
  connect(h2StorageSystem_Compressed.P_comp, systemElectrolyzer.CompressorPower) annotation (Line(
      points={{62.8,-20.8},{70,-20.8},{70,-88},{-7.8,-88},{-7.8,-79.4}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(systemPEMFC.P_FC_tot, add.u1) annotation (Line(
      points={{-12,80.8},{-12,88},{52,88},{52,68},{60,68}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(systemElectrolyzer.P_electrolyzer_tot, add.u2) annotation (Line(
      points={{-6,-38.8},{-6,28},{52,28},{52,56},{60,56}},
      color={0,135,135},
      pattern=LinePattern.Dash));
  connect(add.y, P_HESS) annotation (Line(points={{83,62},{106,62}}, color={0,0,127}));
  connect(P_FC.y, systemPEMFC.P_el_set) annotation (Line(points={{-55,64},{-54,64.8},{-41.6,64.8}}, color={0,0,127}));
  connect(P_EL.y, systemElectrolyzer.P_el_set) annotation (Line(points={{-53,-60},{-40.8,-60}},                     color={0,0,127}));
  connect(systemElectrolyzer.T_environment, T_environment) annotation (Line(points={{-40.6,-43.4},{-80,-43.4},{-80,30},{-110,30}}, color={0,0,127}));
  connect(T_environment, systemPEMFC.T_environment) annotation (Line(points={{-110,30},{-80,30},{-80,76},{-42,76}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={255,213,170},
          fillPattern=FillPattern.Solid), Text(
          extent={{-90,42},{90,-36}},
          textColor={0,0,0},
          textString="HESS
350 bar")}),                                                     Diagram(coordinateSystem(preserveAspectRatio=false)));
end HESS_Compressed;
