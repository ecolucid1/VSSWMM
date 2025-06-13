# VSSWMM
A single python script that reproduces the means and methods of the SWMM engine to predict runoff from a catchment using real-time weather feeds.
# Lightweight SWMM Runoff Predictor

A single-file Python script that predicts real-time stormwater runoff from a single catchment area using the core hydrology calculations from the official EPA Storm Water Management Model (SWMM 5).

This script connects to the OpenWeatherMap API to fetch live precipitation data and applies scientifically validated methods to provide a continuous, real-time estimate of runoff in cubic feet per second (cfs).

## Core Concepts

This tool is a Python translation and simplification of the powerful and complex EPA SWMM 5 C-language engine. It is designed to be lightweight and easy to run for quick, real-time predictions without needing a full SWMM input file or installation.

The script faithfully implements two key SWMM 5 methods:

1.  **Nonlinear Reservoir Model:** Surface runoff is calculated by treating the catchment surface as a nonlinear reservoir where the outflow is governed by the Manning equation.
2.  **Horton Infiltration:** Infiltration into pervious land surfaces is modeled using the widely-accepted Horton method, which simulates the decay of infiltration rate as the soil becomes saturated.
3.  **Adaptive ODE Solver:** To ensure numerical accuracy, the change in surface water depth over time is solved using an adaptive 5th-order Runge-Kutta method, just as in the original SWMM engine.

## Features

* **Real-Time Data:** Fetches live, hourly precipitation data from the OpenWeatherMap "One Call" API.
* **SWMM 5 Calculations:** Utilizes the core hydrology algorithms from the EPA SWMM 5 engine.
* **Configurable Catchment:** Easily configure the physical characteristics of your watershed (area, slope, imperviousness, etc.) within the script.
* **Self-Contained:** All logic is contained in a single Python file with a minimal dependency (`requests`).
* **Educational:** Provides a clear, object-oriented Python implementation of complex hydrological processes, making it a great learning tool.

## How to Use

Follow these steps to get the script running.

### 1. Prerequisites

* Python 3.6+
* The `requests` library. If you don't have it, install it via pip:
    ```bash
    pip install requests
    ```

### 2. Get an API Key

* Sign up for a free account at [OpenWeatherMap](https://openweathermap.org/).
* Navigate to your user dashboard and find the **"API keys"** tab.
* Copy your default API key.

### 3. Configure the Script

* Open the `swmm_runoff_py.py` script in a text editor.
* **Set your API Key:** On line 228, replace `"YOUR_API_KEY"` with the key you copied from OpenWeatherMap.

    ```python
    API_KEY = "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    ```

* **Set Your Location:** On lines 231-232, change the `LATITUDE` and `LONGITUDE` to your area of interest.

    ```python
    LATITUDE = 34.0522    # e.g., Los Angeles
    LONGITUDE = -118.2437
    ```
* **(Optional) Configure the Watershed:** On lines 238-249, you can adjust the `watershed_params` dictionary to match the characteristics of your catchment. The default values represent a typical 20-acre, medium-density residential area.

### 4. Run the Script

* Save your changes to the script.
* Run the script from your terminal:
    ```bash
    python swmm_runoff_py.py
    ```

The script will start fetching data and printing the predicted runoff every 60 seconds.

## Example Output

The terminal output will look like this, updating with each time step:


--- Lightweight SWMM Runoff Predictor --- Simulating runoff for 'ResidentialArea1'... Press Ctrl+C to stop.
Time | Rainfall (in/hr) | Runoff (cfs)
14:30:15   | 0.0000               | 0.0000
14:31:15   | 0.0000               | 0.0000
14:32:16   | 0.1500               | 0.8532
14:33:16   | 0.1500               | 1.2104
...


## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

## Acknowledgments

This script is a Python implementation based on the original **EPA SWMM 5** source code. Credit and thanks go to the original authors and contributors at the U.S. Environmental Protection Agency (EPA) for their extensive work on this important public domain model.

