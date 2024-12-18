within H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.Specifications;
record GinerELX5kW "DESL Giner Electrolyzer 5.5kW system specific parameters"

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

  //Record containing the data of the DESL 5kW Giner Inc. PEMEL stack

  // _____________________________________________
  //
  //          Imports and Class Hierarchy
  // _____________________________________________

  extends H2Microgrid_TransiEnt.ElectrolyzerBoPSystem.Electrolyzer.Specifications.Base5kWElectrolyzerL2Specification(
    n_cells=20,
    alpha_an=0.75818,
    E_exc=47290.321,
    E_pro=14202.4525,
    i_dens_0_an_std=2.1014e-5,
    mem_conductivity_ref=7.3485,
    t_mem=183e-6,
    PEM_area=50e-4,
    i_el_n=150,
    P_el_n=5e3,
    R_th=0.529056,
    C_th=54038.66667/6,
    P_el_max=1*P_el_n,
    P_el_pump = 285.7,
    eta_pumpmotor=0.51,
    V_flow_water=0.00005,
    Delta_p_pump=9.2*100000,
    T_op_max=273.15 + 50.5,
    T_cool_set=273.15 + 50,
    Q_flow_cool_max=1600);
  annotation (Documentation(info="<html>
<h4><span style=\"color: #008000\">1. Purpose of model</span></h4>
<p>Record containing system specifications of a Giner configuration of a Giner 5kW Electrolyzer system</p>
<h4><span style=\"color: #008000\">2. Level of detail, physical effects considered, and physical insight</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">3. Limits of validity </span></h4>
<p>(no remarks) </p>
<h4><span style=\"color: #008000\">4. Interfaces</span></h4>
<p>(no remarks) </p>
<h4><span style=\"color: #008000\">5. Nomenclature</span></h4>
<p>(no remarks) </p>
<h4><span style=\"color: #008000\">6. Governing Equations</span></h4>
<p>(no remarks)</p>
<h4><span style=\"color: #008000\">7. Remarks for Usage</span></h4>
<p>(no remarks) </p>
<h4><span style=\"color: #008000\">8. Validation</span></h4>
<p>(no remarks) </p>
<h4><span style=\"color: #008000\">9. References</span></h4>
<p>Source: modified from TransiEnt library original model</p>
<p><br>Z. Abdin, E. MacA. Gray, and C.J. Webb. Modelling and simulation of a proton exchange membrane (PEM) electrolyzer cell. International Journal of Hydrogen Energy, 40(39):13243-13257, 2015. doi:10.1016/j.ijhydene.2015.07.129. </p>
<p>Manuel Espinosa-L&oacute;pez, Philippe Baucour, Serge Besse, Christophe Darras, Raynal Glises, Philippe Poggi, Andr&eacute; Rakotondrainibe, and Pierre Serre-Combe. Modelling and experimental validation of a 46 kW PEM high pressure water electrolyser. Renewable Energy, 119, pp. 160-173, 2018. doi: 10.1016/J.RENENE.2017.11.081.</p>
<h4><span style=\"color: #008000\">10. Version History</span></h4>
<p>Created by John Webster (jcwebste@edu.uwaterloo.ca) October 2018.</p>
</html>"));
end GinerELX5kW;
