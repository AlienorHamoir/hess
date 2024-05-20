within H2Microgrid_TransiEnt.StorageSystem;
block TankSOC "Block that computes the SOC of a tank based on current pressure"
  parameter Modelica.Units.SI.Pressure maxPressure = 30e5 "Maximum pressure the tank can withstand";

  Modelica.Blocks.Interfaces.RealOutput tankSOC
    annotation (Placement(transformation(extent={{96,-10},{116,10}})));
  Modelica.Blocks.Interfaces.RealInput currentTankPressure
    annotation (Placement(transformation(extent={{-122,-22},{-82,18}})));
equation
  tankSOC = currentTankPressure / maxPressure;
  annotation (Icon(graphics={Text(
          extent={{-52,20},{50,-22}},
          textColor={102,44,145},
          textString="TankSOC"), Rectangle(extent={{-56,18},{56,-18}},
            lineColor={102,44,145})}));
end TankSOC;
