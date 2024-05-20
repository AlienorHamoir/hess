within H2Microgrid_TransiEnt.StorageSystem;
model H2Tank_Pressure "Pressure model of a storage tank (Rullo, 2016)"

  TransiEnt.Basics.Interfaces.General.PressureOut p_tank annotation (Placement(transformation(extent={{60,-20},{98,18}}), iconTransformation(extent={{60,-20},{98,18}})));
  TransiEnt.Basics.Interfaces.General.MassFlowRateIn H2massFlowElectrolyzer annotation (Placement(transformation(extent={{-100,20},{-60,60}}), iconTransformation(extent={{-100,20},{-60,60}})));
  TransiEnt.Basics.Interfaces.General.MassFlowRateIn H2massFlowFC annotation (Placement(transformation(extent={{-100,-60},{-60,-20}}), iconTransformation(extent={{-100,-60},{-60,-20}})));

  parameter Real H2MolecularWeight = 1.00784;
  Real FH2_electrolyzer = H2massFlowElectrolyzer / H2MolecularWeight;
  Real FH2_FC = H2massFlowFC / H2MolecularWeight;

  parameter Real R=8.314;
  parameter Modelica.Units.SI.Temperature T_tank=273.15+10;
  parameter Modelica.Units.SI.Volume V_tank = 0.1;
  replaceable parameter Modelica.Units.SI.Pressure p_start = 17e5;

  Real p_tankPa(start=p_start); // tank pressure in Pa

equation
  der(p_tankPa) = R*T_tank/V_tank*(FH2_electrolyzer-FH2_FC); // tank pressure derivative in Pa zith ideal gas law
  p_tank = 1*p_tankPa; // tank pressure in bar

end H2Tank_Pressure;
