within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.H2DryingSubsystem;
model H2DryingSubsystem "Subsystem in which H2 gets dried and excess water is absorbed by a sink"

  replaceable parameter TILMedia.VLEFluidTypes.TILMedia_SplineWater water;
  replaceable parameter TransiEnt.Basics.Media.Gases.VLE_VDIWA_SG4_var medium_gas;

  Modelica.Units.SI.MassFlowRate m_flow_CH4_in=dryer.gasPortIn.m_flow*inStream(dryer.gasPortIn.xi_outflow[1]);
  Modelica.Units.SI.MassFlowRate m_flow_CH4_out=-(dryer.gasPortOut.m_flow*dryer.gasPortOut.xi_outflow[1]);
  Modelica.Units.SI.MassFlowRate m_flow_CO2_in=dryer.gasPortIn.m_flow*inStream(dryer.gasPortIn.xi_outflow[2]);
  Modelica.Units.SI.MassFlowRate m_flow_CO2_out=-(dryer.gasPortOut.m_flow*dryer.gasPortOut.xi_outflow[2]);
  Modelica.Units.SI.MassFlowRate m_flow_H2O_in=dryer.gasPortIn.m_flow*inStream(dryer.gasPortIn.xi_outflow[3]);
  Modelica.Units.SI.MassFlowRate m_flow_H2O_out=-(dryer.gasPortOut.m_flow*dryer.gasPortOut.xi_outflow[3] + dryer.fluidPortOut.m_flow);
  Modelica.Units.SI.MassFlowRate m_flow_H2_in=dryer.gasPortIn.m_flow*(1 - sum(inStream(dryer.gasPortIn.xi_outflow)));
  Modelica.Units.SI.MassFlowRate m_flow_H2_out=-(dryer.gasPortOut.m_flow*(1 - sum(dryer.gasPortOut.xi_outflow)));
//   Modelica.Units.SI.MassFlowRate m_flow_CO_in=dryer.gasPortIn.m_flow*inStream(dryer.gasPortIn.xi_outflow[5]);
//   Modelica.Units.SI.MassFlowRate m_flow_CO_out=-(dryer.gasPortOut.m_flow*dryer.gasPortOut.xi_outflow[5]);
//   Modelica.Units.SI.MassFlowRate m_flow_N2_in=dryer.gasPortIn.m_flow*(1 - sum(inStream(dryer.gasPortIn.xi_outflow)));
//   Modelica.Units.SI.MassFlowRate m_flow_N2_out=-(dryer.gasPortOut.m_flow*(1 - sum(dryer.gasPortOut.xi_outflow)));

  Dryer_L1 dryer(
    medium_gas=medium_gas,
    medium_water=water,
    pressureLoss=20000) annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

  TransiEnt.Components.Boundaries.FluidFlow.BoundaryVLE_pTxi
                                                           sink_water(
    variable_p=false,
    medium=water,
    boundaryConditions(showData=false))                                                               annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=270,
        origin={0,40})));
  TransiEnt.Basics.Interfaces.Gas.RealGasPortOut gasPortOut(Medium=medium_gas) annotation (Placement(transformation(extent={{90,-10},{110,10}})));
  TransiEnt.Basics.Interfaces.Gas.RealGasPortIn gasPortIn(Medium=medium_gas) annotation (Placement(transformation(extent={{-108,-10},{-88,10}})));

  outer TransiEnt.SimCenter
                  simCenter annotation (Placement(transformation(extent={{-78,54},{-58,74}})));
equation
  connect(sink_water.fluidPortIn, dryer.fluidPortOut) annotation (Line(
      points={{0,30},{0,10}},
      color={175,0,0},
      thickness=0.5));
  connect(dryer.gasPortOut, gasPortOut) annotation (Line(
      points={{10,0},{100,0}},
      color={255,255,0},
      thickness=1.5));
  connect(dryer.gasPortIn, gasPortIn) annotation (Line(
      points={{-10,0},{-98,0}},
      color={255,255,0},
      thickness=1.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end H2DryingSubsystem;
