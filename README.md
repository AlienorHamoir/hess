# README

### Summary

This repository contains the Modelica models for a hybrid microgrid and a Model Predictive Control (MPC) algorithm, developed in Python, which optimized the microgrid operations based on the total cast. It is designed to run in co-simulation with a model exported as an FMU. The hybrid microgrid consists of a Hydrogen Energy Storage System, a battery, PV generation, and a load.
The MCP works with forecasts for electricity prices, PV generation and load as inputs. Different MPC versions have been developped, encompassing different modelling levels of detail, in order to assess the impact of the modelling accuracy on the optimizer performances.

**Version**

  - Modelica: 4.0X
  - Python: 3.9

### Setup

#### Modelica Models

The Modelica models (`.mos` and `.mo` files) must run on Dymola. These models are developed using components from the TransiEnt, SCooDER, and Buildings libraries. 

- **TransiEnt Library**: Install by following instructions on [TransiEnt GitHub](https://github.com/TransiEnt-official/transient-lib).

- **Loading TransiEnt Library in Dymola**: Use `loadTransiEnt.mos` from this repository, adjusting the file paths as necessary. Ensure compatibility by using the modified `simCenter` file found in `H2Microgrid_TransiEnt/Resources/simC_modif.txt`.

The `H2Microgrid_TransiEnt` folder contains `.mo` and `.mos` files for different parts, subparts, and submodels of the hybrid microgrid, including the electrolyzer, fuel cell and hydrogen tank subsystems, as well as the assembled hdyrogen energy storage systems and the microgrid.

#### MPC Algorithm and Co-simulation Setup

- The Jupyter Notebook for co-simulation with the MPC algorithm can be found at `MPC-FMU/Co-simulation_MPC.ipynb`.

- `MPC-FMU/MPC.ipynb` serves as a test notebook for various MPC versions.

- **Forecasts**: `MPC-FMU/forecasts` contains saved `.csv` files with scenarios for solar production, load, and electricity prices (e.g., winter/summer, weekdays/weekends). To visualize forecasts, use `forecasts_view.ipynb`.

#### FMU and Co-simulation Environment

- The latest FMU export, `H2Microgrid_0TransiEnt_HybridMicrogrid_H2Microgrid_0HP.fmu`, is available in `MPC-FMU/FMU`.

- For the MPC to operate in co-simulation, a functional mock-up unit (FMU) of the hybrid microgrid Dymola model is used. The FMU must be downloaded to the local machine.

- **FMI-MLC Package**: Required to run the MPC with the FMU in co-simulation. This package requires a specific environment accessible through `StartJupyter_generic.bat`.
  
  - **Docker**: Docker Desktop must be installed locally to launch the necessary environment container.

#### Optimization: Gurobi and CasADi

- The MPC relies on the Gurobi solver for optimization. A local Gurobi license must be installed to use this solver. CasADi is used as the optimization toolbox.

- **Important**: The Gurobi solver is not currently configured on the Docker containers Linux environment, so full MPC functionality may be unavailable.

#### Recommended Setup and Configuration

- **IDE**: Running Python code in VSCode is recommended to avoid kernel crashes. Avoid running code in Jupyter Notebooks or Google Colab.

- **Docker Configuration**: Ensure the Docker container is connected to the Linux machine hosting Docker.

### Usage

1. **Open the Docker environment** by executing the StartJupyter_generic.bat file or by running the following command in your terminal:
         ./StartJupyter_generic.bat

2. **Open the Jupyter Notebook for Co-simulation**:
   
   - Navigate to `MPC-FMU/Co-simulation_MPC.ipynb`.
   
   - Follow the instructions and run the cells in order to execute the MPC algorithm.

3. **View Results**:
   
   - Output files, including optimal operational profiles, will be saved in the `results` directory.
   
   - Check the generated `.pdf` plots and `.csv` files for analysis.


#### Connecting to the Docker Container from VSCode

To work with the Docker container directly from Visual Studio Code, follow these steps:

1. **Install Required VSCode Extensions**:
   
   - **Remote - Containers**: This extension allows you to connect to and develop inside a Docker container directly from VSCode. Install it from the [VSCode Extensions Marketplace](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers).

2. **Connect to the Container**:
   
   - **Start Docker Desktop** and ensure the container you want to connect to is running.
   
   - In VSCode, open the **Command Palette** (press `Cmd+Shift+P` on macOS or `Ctrl+Shift+P` on Windows/Linux).
   
   - Type and select **Remote-Containers: Attach to Running Container...**.
   
   - Choose your container from the list of running containers.

3. **Open and Edit Files in the Container**:
   
   - Once connected, VSCode will open a new window where you can browse, edit, and execute code within the container environment.
   
   - You can access all files and folders available inside the container, allowing seamless integration with the MPC and FMU models.

4. **Installing Additional Dependencies in the Container**:
   
   - If you need to install any additional Python packages or dependencies, you can open a terminal inside VSCode (from **Terminal > New Terminal**) and run commands directly within the container environment.

**Note**: Ensure Docker Desktop has sufficient permissions to access the directories and files you’ll be working with, particularly on macOS where you may need to enable directory access under **Docker > Settings > Resources > File Sharing**.

This setup allows you to develop and debug within the containerized environment without switching between the host and container, providing a more integrated workflow.

### Workflow

0. **Parameters and simualtion features configurations:** The time horizon and the timestep of the simulation can be adjusted in the code before the forecast extraction from the FMU. The start and end date and time of the simualtion can be chosen as well.

1. **Data Extraction**: The FMU is initially read to extract PV production and load data, which are then stored and used in the MPC.

2. **MPC Versions**: Different MPC versions with varying levels of modeling detail are available, each defined in separate functions. One version of the MPC should be chosen for each run, depending on the test and results. The desired version to be called shall be changed in the code. Details on the different MPC versions developed can be found in 'Alienor Hamoir Master Thesis - Project Documentation'.

3. **Results**: 

   - Outputs include the optimal predicted operational profile for the hydrogen energy storage system and battery, grid exchange data, and operational cost.
   
   - The results are saved as `.pdf` files for plots and `.csv` files (to be added) for raw data.

To run tests, execute each cell in order.

### Results

The simulation generates:

- **PDF Plots**: Visualization of operational profiles.

- **CSV Files**: Raw data for further analysis, including operational profiles for the hydrogen energy storage system, battery, and grid exchange data. (to be added)

### Folder Structure

- `H2Microgrid_TransiEnt/`: Contains Modelica files and submodels for the microgrid.

- `MPC-FMU/`: Holds the Jupyter Notebooks for MPC simulation and forecast data.

- `resources/`: Contains weather, load, and FMU files required for the simulations.

**Note**: The `simulationresults` folder can be ignored.

### Resources

The `resources` folder contains:

- **Weather Files**: Files for various U.S. zones, used by the Dymola model.

- **Load Files**: Includes `.idf` files and formatted `.txt` files. The `commercial_smallOffice_LA.txt` file is used in analysis.

- **FMU Files**: Latest FMU export is available in `fmus`.


### Dependencies

- Modelica 4.0
- Python 3.9
- Dymola
- TransiEnt, SCooDER, and Buildings libraries for Modelica
- Gurobi and CasADi (Python packages for optimization)
- FMI-MLC Package (requires Docker)

### Additional References

- [TransiEnt Documentation](https://github.com/TransiEnt-official/transient-lib)
- [CasADi Documentation](https://web.casadi.org/)
- [Dymola Documentation](https://www.3ds.com/products-services/catia/products/dymola/)
- [Alienor Hamoir Master Thesis - Project Documentation](https://drive.google.com/file/d/1P-CO3N2X6Iqkec_dlwpKemOehupPUwS_/view?usp=share_link)


### Contribution Guidelines

- **Code Reviews**: Please follow standard practices for code quality and readability.
- **Testing**: Test each major change to ensure compatibility with the existing setup.
- **Documentation**: Ensure all new contributions are well-documented.
- **Other guidelines**


### Database configuration

### Deployment instructions


### Developed By

Alienor Hamoir  
Contact: [alienor.hamoir@hotmail.com](mailto:alienor.hamoir@hotmail.com)
