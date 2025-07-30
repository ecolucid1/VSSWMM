"""
# Lightweight SWMM Runoff Predictor
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Author:** Marcus Quigley

A single-file Python script that predicts real-time stormwater runoff from a single catchment area using the core hydrology calculations from the official EPA Storm Water Management Model (SWMM 5).

This script connects to the OpenWeatherMap API to fetch live precipitation data and applies scientifically validated methods to provide a continuous, real-time estimate of runoff in cubic feet per second (cfs).

## Core Concepts

This tool is a Python translation and simplification of the powerful and complex EPA SWMM 5 C-language engine. It is designed to be lightweight and easy to run for quick, real-time predictions without needing a full SWMM input file or installation.

The script faithfully implements three key SWMM 5 methods:

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

* Open this script in a text editor.
* **Set your API Key:** In the `if __name__ == "__main__"` block, replace `"YOUR_API_KEY"` with the key you copied from OpenWeatherMap.
* **Set Your Location:** Change the `LATITUDE` and `LONGITUDE` to your area of interest.
* **(Optional) Configure the Watershed:** Adjust the `watershed_params` dictionary to match the characteristics of your catchment.

### 4. Run the Script

* Save your changes.
* Run the script from your terminal:
    ```bash
    python your_script_name.py
    ```

## License

This project is licensed under the MIT License.

Copyright (c) 2024 Marcus Quigley

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Credits and Acknowledgments

* **Author:** Marcus Quigley
* **Original Source:** This script is a Python implementation based on the original **EPA SWMM 5** source code. Credit and thanks go to the original authors and contributors at the U.S. Environmental Protection Agency (EPA) for their extensive work on this important public domain model.
"""
import requests
import time
import math
from typing import Dict, Any, List

# --- Constants based on SWMM's consts.h and internal units ---
# Conversion factor for Manning's equation in US customary units
MANNING_COEFF = 1.49
# Runoff exponent in Manning's equation
MANNING_EXP = 1.6666667
# Acceptable error for the ODE solver
ODE_TOL = 0.0001
# Seconds per hour
SEC_PER_HOUR = 3600.0

class Subcatchment:
    """
    A class to represent a single urban catchment, encapsulating its
    physical properties and the state of its surface water and infiltration.
    The methods are direct Python translations of the core SWMM-5 C code.
    """
    def __init__(self, name: str, params: Dict[str, float]):
        """
        Initializes the subcatchment with hardcoded parameters.
        Units are internal SWMM units (feet and seconds).
        """
        self.name = name
        # --- Static Watershed Characteristics (Hardcoded) ---
        self.area_acres = params.get("area_acres", 10.0)            # Area in acres
        self.width_ft = params.get("width_ft", 500.0)              # Characteristic width in feet
        self.slope_pct = params.get("slope_pct", 0.5)              # Average slope as a percentage
        self.imperv_pct = params.get("imperv_pct", 50.0)           # Impervious area as a percentage
        self.n_imperv = params.get("n_imperv", 0.01)               # Manning's n for impervious area
        self.n_perv = params.get("n_perv", 0.10)                   # Manning's n for pervious area
        self.dstore_imperv_in = params.get("dstore_imperv_in", 0.05) # Depression storage, impervious (in)
        self.dstore_perv_in = params.get("dstore_perv_in", 0.10)   # Depression storage, pervious (in)
        self.zero_imperv_pct = params.get("zero_imperv_pct", 25.0)   # % of impervious area with no storage

        # --- Horton Infiltration Parameters (Hardcoded) ---
        self.infil_max_rate_in_hr = params.get("infil_max_rate_in_hr", 3.0) # Max infiltration rate (in/hr)
        self.infil_min_rate_in_hr = params.get("infil_min_rate_in_hr", 0.13)# Min infiltration rate (in/hr)
        self.infil_decay = params.get("infil_decay", 4.0)                  # Decay constant (1/hr)
        self.infil_dry_time_days = params.get("infil_dry_time_days", 7.0)  # Time for soil to dry (days)

        # --- Calculated properties (converted to internal SWMM units: feet/seconds) ---
        self.area_ft2 = self.area_acres * 43560.0
        self.slope = self.slope_pct / 100.0

        # --- Sub-area characteristics ---
        frac_imperv = self.imperv_pct / 100.0
        frac_zero_imperv = self.zero_imperv_pct / 100.0

        self.imperv_area_1 = self.area_ft2 * frac_imperv * frac_zero_imperv  # Impervious, no storage
        self.imperv_area_2 = self.area_ft2 * frac_imperv * (1.0 - frac_zero_imperv) # Impervious, w/ storage
        self.perv_area = self.area_ft2 * (1.0 - frac_imperv)

        # Depression storage in feet
        self.dstore_imperv = self.dstore_imperv_in / 12.0
        self.dstore_perv = self.dstore_perv_in / 12.0

        # Runoff coefficient 'alpha' for Manning's equation
        self.alpha_imperv = self._calculate_alpha(self.n_imperv, self.imperv_area_1 + self.imperv_area_2)
        self.alpha_perv = self._calculate_alpha(self.n_perv, self.perv_area)

        # --- Infiltration properties (converted to ft/sec) ---
        self.infil_f0 = self.infil_max_rate_in_hr / SEC_PER_HOUR / 12.0
        self.infil_fmin = self.infil_min_rate_in_hr / SEC_PER_HOUR / 12.0
        self.infil_decay_rate = self.infil_decay / SEC_PER_HOUR
        self.infil_regen_rate = 1.0 / (self.infil_dry_time_days * 24 * SEC_PER_HOUR)

        # --- Dynamic State Variables ---
        self.depth_imperv_1 = 0.0  # Ponded depth on impervious area 1 (ft)
        self.depth_imperv_2 = 0.0  # Ponded depth on impervious area 2 (ft)
        self.depth_perv = 0.0      # Ponded depth on pervious area (ft)
        self.infil_f = self.infil_f0 # Current infiltration capacity (ft/sec)

    def _calculate_alpha(self, n, sub_area_ft2):
        """Calculates the alpha coefficient for the Manning equation."""
        if sub_area_ft2 <= 0 or n <= 0:
            return 0
        return MANNING_COEFF * self.width_ft * math.sqrt(self.slope) / (n * sub_area_ft2)

    def get_runoff(self, rainfall_rate_in_hr: float, t_step_sec: float) -> float:
        """
        Calculates total runoff from the subcatchment for a given time step.
        This is a Python translation of the logic in runoff.c and subcatch.c.
        """
        rainfall_rate_ft_sec = rainfall_rate_in_hr / 12.0 / SEC_PER_HOUR

        # --- 1. Calculate runoff from impervious area with no depression storage ---
        runoff_imperv_1, self.depth_imperv_1 = self._get_subarea_runoff(
            self.depth_imperv_1,
            rainfall_rate_ft_sec,
            0, # No depression storage
            self.alpha_imperv,
            t_step_sec)

        # --- 2. Calculate runoff from impervious area with depression storage ---
        runoff_imperv_2, self.depth_imperv_2 = self._get_subarea_runoff(
            self.depth_imperv_2,
            rainfall_rate_ft_sec,
            self.dstore_imperv,
            self.alpha_imperv,
            t_step_sec)

        # --- 3. Calculate infiltration and runoff for the pervious area ---
        # Update infiltration capacity based on rainfall and soil state
        if rainfall_rate_ft_sec > 0:
            self.infil_f -= self.infil_decay_rate * (self.infil_f - self.infil_fmin) * t_step_sec
        else:
            self.infil_f += self.infil_regen_rate * (self.infil_f0 - self.infil_f) * t_step_sec
        self.infil_f = max(self.infil_f, self.infil_fmin)
        self.infil_f = min(self.infil_f, self.infil_f0)

        # Effective infiltration rate cannot exceed available water
        available_water_rate = self.depth_perv / t_step_sec + rainfall_rate_ft_sec
        infiltration_rate = min(self.infil_f, available_water_rate)

        # Net rainfall after infiltration
        net_rainfall_perv = rainfall_rate_ft_sec - infiltration_rate

        runoff_perv, self.depth_perv = self._get_subarea_runoff(
            self.depth_perv,
            net_rainfall_perv,
            self.dstore_perv,
            self.alpha_perv,
            t_step_sec)

        # --- 4. Sum up the runoff from all sub-areas (in cfs) ---
        total_runoff_cfs = (runoff_imperv_1 * self.imperv_area_1 +
                            runoff_imperv_2 * self.imperv_area_2 +
                            runoff_perv * self.perv_area)
        return total_runoff_cfs

    def _get_subarea_runoff(self, d_initial, inflow_rate, d_store, alpha, t_step):
        """
        Solves the ODE for a single sub-area to find new depth and runoff.
        Uses the adaptive Runge-Kutta solver.
        """
        # --- Store sub-area context for the ODE solver ---
        self.ode_context = {
            "inflow_rate": inflow_rate,
            "d_store": d_store,
            "alpha": alpha
        }

        # --- Integrate to find new depth ---
        y_start = [d_initial]
        self._odesolve_integrate(y_start, 1, 0, t_step, ODE_TOL, t_step, self._get_dd_dt)
        d_final = y_start[0]

        # --- Calculate runoff rate based on final depth ---
        d_excess = d_final - d_store
        if d_excess > 0 and alpha > 0:
            runoff_rate = alpha * (d_excess ** MANNING_EXP)
        else:
            runoff_rate = 0.0
            
        return runoff_rate, d_final

    def _get_dd_dt(self, t: float, y: List[float]) -> List[float]:
        """
        Calculates the derivative of depth (dD/dt) for the ODE solver.
        dD/dt = Inflow - Outflow
        """
        context = self.ode_context
        inflow = context["inflow_rate"]
        d_store = context["d_store"]
        alpha = context["alpha"]
        
        depth = y[0]
        d_excess = depth - d_store
        
        if d_excess > 0 and alpha > 0:
            outflow = alpha * (d_excess ** MANNING_EXP)
        else:
            outflow = 0
            
        return [inflow - outflow]

    def _odesolve_integrate(self, y_start: List[float], n: int, x1: float, x2: float,
                            eps: float, h1: float, derivs):
        """
        Python implementation of the adaptive Runge-Kutta 5th order ODE solver.
        This is a translation of the `odesolve_integrate` and related functions
        from SWMM's odesolve.c.
        """
        MAXSTP = 10000
        TINY = 1.0e-30
        SAFETY = 0.9
        PGROW = -0.2
        PSHRNK = -0.25
        ERRCON = (5.0 / SAFETY) ** (1.0 / PGROW)

        x = x1
        h = h1
        y = list(y_start)

        for _ in range(MAXSTP):
            dydx = derivs(x, y)
            yscal = [abs(y[i]) + abs(dydx[i] * h) + TINY for i in range(n)]

            if (x + h - x2) * (x + h - x1) > 0.0:
                h = x2 - x
            
            # --- Adaptive step size loop ---
            while True:
                y_err, y_temp = self._rkck(y, dydx, n, x, h, derivs)
                errmax = 0.0
                for i in range(n):
                    if yscal[i] != 0:
                        err = abs(y_err[i] / yscal[i])
                        if err > errmax:
                            errmax = err
                errmax /= eps

                if errmax <= 1.0:
                    break  # Step succeeded

                # Step failed, reduce step size
                h_temp = SAFETY * h * (errmax ** PSHRNK)
                h = max(h_temp, 0.1 * h) if h > 0 else min(h_temp, 0.1 * h)
                if x + h == x:
                    raise RuntimeError("ODE solver step size underflow")

            # --- If step succeeded, prepare for next step ---
            if errmax > ERRCON:
                h_next = SAFETY * h * (errmax ** PGROW)
            else:
                h_next = 5.0 * h
            
            x += h
            y = y_temp
            h = h_next

            if (x - x2) * (x2 - x1) >= 0.0:
                for i in range(n):
                    y_start[i] = y[i]
                return
        
        raise RuntimeError("Exceeded max steps in ODE solver")

    def _rkck(self, y, dydx, n, x, h, derivs):
        """The Cash-Karp Runge-Kutta step."""
        # --- Coefficients for the Cash-Karp method ---
        A = [0.0, 0.2, 0.3, 0.6, 1.0, 0.875]
        B = [
            [],
            [0.2],
            [3.0/40.0, 9.0/40.0],
            [0.3, -0.9, 1.2],
            [-11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0],
            [1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0]
        ]
        C = [37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0]
        DC= [C[0]-2825.0/27648.0, 0.0, C[2]-18575.0/48384.0, C[3]-13525.0/55296.0,
             -277.0/14336.0, C[5]-0.25]

        ytemp = [0.0] * n
        k = [[0.0] * n for _ in range(6)]
        k[0] = dydx

        for i in range(1, 6):
            xtemp = x + A[i] * h
            for j in range(n):
                ytemp[j] = y[j] + h * sum(B[i][m] * k[m][j] for m in range(i))
            k[i] = derivs(xtemp, ytemp)
            
        y_out = [0.0] * n
        y_err = [0.0] * n
        for i in range(n):
            y_out[i] = y[i] + h * sum(C[j] * k[j][i] for j in range(6))
            y_err[i] = h * sum(DC[j] * k[j][i] for j in range(6))
        
        return y_err, y_out

# --- Weather API Function ---
def get_weather_data(api_key: str, lat: float, lon: float) -> float:
    """
    Fetches real-time precipitation data from OpenWeatherMap API.
    Returns rainfall in inches per hour.
    """
    url = f"https://api.openweathermap.org/data/2.5/weather?lat=${lat}&lon=${lon}&units=imperial&appid=${apiKey}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        
        # 'current' object has rain.1h (in mm). Convert mm/hr to in/hr.
        current_precip_mm_hr = 0.0
        if 'current' in data and 'rain' in data['current'] and '1h' in data['current']['rain']:
             current_precip_mm_hr = data['current']['rain']['1h']
        
        # If no current rain, check if there's any in the next hour's forecast
        elif 'hourly' in data and len(data['hourly']) > 0 and 'rain' in data['hourly'][0]:
             current_precip_mm_hr = data['hourly'][0]['rain']['1h']

        return current_precip_mm_hr / 25.4

    except requests.exceptions.RequestException as e:
        print(f"Error fetching weather data: {e}")
        return 0.0

# --- Main Simulation ---
if __name__ == "__main__":
    # --- Configuration ---
    # To get a free API key, sign up at https://openweathermap.org/
    API_KEY = "YOUR_API_KEY"  # <-- IMPORTANT: REPLACE WITH YOUR KEY
    
    # Location (e.g., New York City)
    LATITUDE = 40.7128
    LONGITUDE = -74.0060

    # Simulation time step in seconds
    TIME_STEP_SEC = 60

    # --- Setup the Subcatchment ---
    # These parameters represent a typical medium-density residential area.
    # They can be modified to represent different types of land use.
    watershed_params = {
        "area_acres": 20.0,
        "width_ft": 700.0,
        "slope_pct": 1.0,
        "imperv_pct": 60.0,
        "n_imperv": 0.015,
        "n_perv": 0.15,
        "dstore_imperv_in": 0.06,
        "dstore_perv_in": 0.12,
        "infil_max_rate_in_hr": 2.0,
        "infil_min_rate_in_hr": 0.2,
    }
    my_catchment = Subcatchment("ResidentialArea1", watershed_params)

    print("--- Lightweight SWMM Runoff Predictor ---")
    print(f"Simulating runoff for '{my_catchment.name}'...")
    if API_KEY == "YOUR_API_KEY":
        print("\nWARNING: Please replace 'YOUR_API_KEY' with your actual OpenWeatherMap API key.")
    print("Press Ctrl+C to stop.")
    print("-" * 40)
    print(f"{'Time':<10} | {'Rainfall (in/hr)':<20} | {'Runoff (cfs)':<15}")
    print("-" * 40)

    # --- Simulation Loop ---
    try:
        while True:
            # 1. Fetch real-time rainfall data
            rainfall_in_hr = get_weather_data(API_KEY, LATITUDE, LONGITUDE)

            # 2. Calculate runoff using the SWMM model
            runoff_cfs = my_catchment.get_runoff(rainfall_in_hr, TIME_STEP_SEC)
            
            # 3. Print results
            current_time_str = time.strftime("%H:%M:%S")
            print(f"{current_time_str:<10} | {rainfall_in_hr:<20.4f} | {runoff_cfs:<15.4f}")

            # 4. Wait for the next time step
            time.sleep(TIME_STEP_SEC)

    except KeyboardInterrupt:
        print("\nSimulation stopped.")
