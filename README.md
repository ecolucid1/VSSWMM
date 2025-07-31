# Lightweight SWMM Runoff Predictor

This project provides two ways to predict real-time stormwater runoff from a single catchment area using the core hydrology calculations from the official EPA Storm Water Management Model (SWMM 5):

1.  **Python Script (`swmm_runoff_py.py`):** A single, self-contained Python script for command-line execution.
2.  **Web Application (`index.html`):** A user-friendly, interactive web interface that runs the same simulation engine directly in your browser.

Both tools connect to the OpenWeatherMap API to fetch live precipitation data and apply scientifically validated methods to provide a continuous, real-time estimate of runoff in cubic feet per second (cfs).

## Core Concepts

This tool is a translation and simplification of the powerful and complex EPA SWMM 5 C-language engine. It is designed to be lightweight and easy to run for quick, real-time predictions without needing a full SWMM input file or installation.

The project faithfully implements three key SWMM 5 methods:

* **Nonlinear Reservoir Model:** Surface runoff is calculated by treating the catchment surface as a nonlinear reservoir where the outflow is governed by the Manning equation.
* **Horton Infiltration:** Infiltration into pervious land surfaces is modeled using the widely-accepted Horton method, which simulates the decay of infiltration rate as the soil becomes saturated.
* **Adaptive ODE Solver:** To ensure numerical accuracy, the change in surface water depth over time is solved using an adaptive 5th-order Runge-Kutta method, just as in the original SWMM engine.

---

## Web Application Interface (`VSSWMM_v2.html`)

The interactive web application provides a visual and intuitive way to run simulations without any setup.

### Web App Features

* **Interactive Controls:** Adjust all watershed and infiltration parameters using intuitive sliders.
* **Clickable Map:** Select your simulation location visually by clicking anywhere on a world map.
* **Live Data Visualization:** A real-time chart plots both rainfall and the calculated runoff, offering immediate visual feedback.
* **Manual Rainfall Mode:** A dedicated "Manual Mode" with a slider allows you to override live weather data to test specific storm scenarios and see how the watershed responds.
* **No Installation Required:** The entire application is self-contained in a single HTML file and runs directly in any modern web browser.

### How to Use the Web App

1.  Open the `index.html` file in a web browser (like Chrome, Firefox, or Edge).
2.  **To use live weather data:**
    * Sign up for a free account at [OpenWeatherMap](https://openweathermap.org/) to get an API key.
    * Enter your API key into the designated field.
    * Use the interactive map to click on your location of interest.
3.  **To simulate a storm:**
    * Click the **"Manual Rainfall"** toggle switch.
    * Use the **"Intensity (in/hr)"** slider to set a desired rainfall rate.
4.  Adjust any of the watershed or infiltration parameters using the sliders in the collapsible sections.
5.  Click **"Start Simulation"**. The results will appear in the chart and table.

---

## Python Script (`swmm_runoff_py.py`)

The original command-line script provides a lightweight way to run the simulation from your terminal.

### Python Script Features

* **Real-Time Data:** Fetches live, hourly precipitation data from the OpenWeatherMap API.
* **SWMM 5 Calculations:** Utilizes the core hydrology algorithms from the EPA SWMM 5 engine.
* **Configurable Catchment:** Easily configure the physical characteristics of your watershed (area, slope, imperviousness, etc.) directly within the script.
* **Self-Contained:** All logic is contained in a single Python file with a minimal dependency (`requests`).
* **Educational:** Provides a clear, object-oriented Python implementation of complex hydrological processes.

### How to Use the Python Script

#### 1. Prerequisites

* Python 3.6+
* The `requests` library. If you don't have it, install it via pip:
    ```bash
    pip install requests
    ```

#### 2. Get an API Key

* Sign up for a free account at [OpenWeatherMap](https://openweathermap.org/).
* Navigate to your user dashboard and find the **"API keys"** tab.
* Copy your default API key.

#### 3. Configure the Script

* Open the `swmm_runoff_py.py` script in a text editor.
* **Set your API Key:** Replace `"YOUR_API_KEY"` with the key you copied.
* **Set Your Location:** Change the `LATITUDE` and `LONGITUDE` to your area of interest.
* **(Optional) Configure the Watershed:** Adjust the `watershed_params` dictionary to match the characteristics of your catchment.

#### 4. Run the Script

* Save your changes and run the script from your terminal:
    ```bash
    python swmm_runoff_py.py
    ```
* The script will start fetching data and printing the predicted runoff to the console.

## License

This project is licensed under the MIT License.

## Acknowledgments

This project is a Python and JavaScript implementation based on the original **EPA SWMM 5** source code. Credit and thanks go to the original authors and contributors at the U.S. Environmental Protection Agency (EPA) for their extensive work on this important public domain model.
