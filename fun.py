"""
Module contains functions to read and calculate discharge (inflow/outflow) and volume concentrations for the generation
of the input boundary condition file for SSIIMM2

"""
from config import *


# Time functions
def modify_time_interval(df_time, flow_array, interval):
    """
    Function re-samples the initial time intervals to a new interval and averages the inflow data, for each watershed,
    for the new time interval.

    It first generates a date list, beginning from the initial time interval to the last, with the desired interval.
    It then does 2 loops, the 1st one is through each new time interval and the 2nd through the original time data
    frame. The original dates are read in order, and therefore the 2nd loop goes through all intervals within the new
    time interval, and averages the flow data.

    NOT USED: Can be used as a guide if something other than averaging for the given time period is to be done

    Args:
    ----------------------------
    :param df_time: data frame, with original time intervals, in datetime format
    :param flow_array: np.array, with inflow data
    :param interval: 1, if re-sample is to daily, 2 if re-sample is to monthly

    :return: data frame, with new time intervals, and np.array, with the flow data averaged for the new time interval.

    Note: Function currently only re-samples to daily and monthly data, and averages the flow data for the new interval.
    If the user wants a different time interval or a different way to manage the flow data, they must do the
    corresponding changes.
    """

    start_date = df_time[0]
    end_date = df_time.iloc[-1]

    if interval == 1:  # Daily (frequency is 'D')
        resampled_dates = pd.to_datetime(pd.date_range(start_date, end_date, freq='D').strftime("%Y%m%d").tolist(),
                                   format="%Y%m%d")
        resampled_flow = np.full((len(resampled_dates), flow_array.shape[1]), 0.0)
        initial = 0
        for d, date in enumerate(resampled_dates):            # Loop through each day to consider
            for i in range(initial, df_time.shape[0]):  # Loop through the original flow/date data
                if date.day != df_time[i].day:
                    resampled_flow[d, :] = np.mean(flow_array[initial:i, :], axis=0)

                    initial = i
                    break
    elif interval == 2:  # monthly (frequency is "MS")
        resampled_dates = pd.to_datetime(pd.date_range(start_date, end_date, freq='MS').strftime("%Y%m%d").tolist(),
                                   format="%Y%m%d")
        resampled_flow = np.full((len(resampled_dates), flow_array.shape[1]), 0.0)
        initial = 0
        for d, date in enumerate(resampled_dates):  # Loop through each day to consider
            for i in range(initial, df_time.shape[0]):  # Loop through the original flow/date data
                if date.month != df_time[i].month:
                    resampled_flow[d, :] = np.mean(flow_array[initial:i, :], axis=0)

                    initial = i
                    break
    else:
        print(f"No resampling was done, since {interval} is not a valid entry")

        return df_time, flow_array

    return resampled_dates, resampled_flow


def resample_time(df_time, flow_array, interval):
    """
    Function re-samples the initial time intervals to a new interval and averages the inflow data, for each watershed,
    for the new time interval.

    It first converts the flow array to a data frame, with the original time data as index. It then uses the df.resample
    tool to resample the data to a given frequency, by averaging the values.

    Args:
    ----------------------------
    :param df_time: data frame Datetime Index or time Series, with original time intervals
    :param flow_array: np.array, with inflow data
    :param interval: 1, if re-sample is to daily, 2 if re-sample is to monthly

    :return: DateTime index, with new time intervals, and np.array, with the flow data averaged for the new
    time interval.

    Note: Function currently only re-samples to daily and monthly data, and averages the flow data for the new interval.
    If the user wants a different frequency, include options for the wanted frequency
    (check: https://regenerativetoday.com/pandas-data_range-function/)
    """

    df_total = pd.DataFrame(data=flow_array, index=df_time)

    if interval == 1:    # Day
        freq = "D"
    elif interval == 2:  # Month
        freq = "MS"
    else:
        print(f"No resampling was done, since {interval} is not a valid entry")
        return df_time, flow_array

    df_resampled = df_total.resample(freq).mean()
    resampled_flow = np.array(df_resampled)
    resampled_dates = pd.to_datetime(df_resampled.index)

    return resampled_dates, resampled_flow


def time_to_seconds(df_time):
    """
    Function calculates the seconds between the beginning of the simulation and each time interval and then transforms
    the differences to seconds.

    Args:
    ------------------------------------
    :param df_time: DataFrame Series or DatetimeIndex, time series with each entry being a time interval for the
    simulation.

    :return: df, time series in seconds.

    Note: if df_time is a DatetimeIndex type, it is transformed to Series type.
    """
    if isinstance(df_time, pd.DatetimeIndex):  # if Datetime index, change to Series
        df_time = df_time.to_series()
    seconds = (df_time - df_time[0]).dt.total_seconds()
    return seconds


# Discharge functions
def extract_discharge(input_df):
    """
    Function reads the raw .b16 discharge output file, in data frame format, and extracts the discharges for the Devoll,
    Holta, Zalli and Skebices sub-catchments. In the input file, columns with names {2, 3, 4, 5} (or indexes
    {6, 7, 8, 9}) corresponds to the catchments {Zalli, Skebices, Holta, Devoll}. The output file needs the following
    sub-catchment order: {Devoll, Holta, Zalli and Skebices}

    Args:
    ----------------------------------------------
    :param input_df: data frame, with the original input data from the .b16 file

    :return: np array, with the discharges for the [Devoll, Holta, Zalli and Skebices] sub catchments, in that order
    """
    output_q = np.full((input_df.shape[0], 4), 0.0)

    output_q[:, 0] = input_df['5']   # Devoll
    output_q[:, 1] = input_df['4']   # Holta
    output_q[:, 2] = input_df['2']   # Zalli
    output_q[:, 3] = input_df['3']   # Skebices

    return output_q


def calculate_outflows_constant_wl(inflow_array, t_capacity):
    """
    Function calculates the outflows based on a constant water level, meaning that inflows = outflows. The outflows are
    through a turbine, with a maximum capacity, and a spillway, in case the inflow exceeds the turbine capacity. The
    total inflow is calculated as the sum of the inflows from each sub-catchment (row-wise sum), and the outflow is
    calculated for each time interval (each row).

    Args:
    --------------------------------------------
    :param inflow_array: np.array, with the inflow data for each sub-catchment (in this case 4). Each column is a
    different sub-catchment and each row corresponds to a different time interval.
    :param t_capacity: float, turbine capacity

    :return: np.array with [inflows (4), turbine outflow, spillway outflow], in that order.
    """
    total_inflow = np.sum(inflow_array, axis=1)

    turbine_outflow = np.where(total_inflow > t_capacity, t_capacity, total_inflow)
    spillway_outflow = np.where(total_inflow > t_capacity, total_inflow - t_capacity, 0)

    # Total flow data
    total_flow = np.c_[inflow_array, turbine_outflow, spillway_outflow]

    return total_flow


def calculate_outflows_seasonal_wl(inflow_array, df_time, sc_path):
    """
    This function is an example of how one could calculate outflow based on water level in the reservoir, using the
    Runge-Kutta method to do the hydrological routing. It needs an initial water level, as well as an equation or a
    discharge curve to associate outflow with water height. It calculates the new water level in each time step and,
    with this water level, it calculates the outflow for the next time step (t+1) based on the data from the current
    time step (t).
    NOTE: For this example, the Q vs h relationship is done using the orifice equation, using a randomly
    chosen diameter and maximum water level, calculated with a function 'alculate_outflow_rk' which does not necessarily
    fit the case. This MUST be modified to fit the real situation and the results should be double-checked, since this
    example is for demonstration purposes only.

    Args:
    ---------------------------------------------
    :param inflow_array: np.array, with the inflow data for each sub-catchment (in this case 4). Each column is a
    different sub-catchment and each row corresponds to a different time interval.
    :param df_time: Data Frame (time Series or DatetimeIndex), with time intervals for the final simulation
    :param sc_path: path where the storage curve is located (in a .txt, tab-delimited file with 2 columns)

    :return:np.array with [inflows (4), turbine outflow, spillway outflow], in that order.
    """
    # save outflow data:
    outflow_array = np.full((inflow_array.shape[0], 2), 0.0)
    h_array = np.full((inflow_array.shape[0], 1), 0.0)

    # time step:
    dt = (df_time[1] - df_time[0]).total_seconds()

    # Total inflow: sum of all inflows
    total_inflow = np.sum(inflow_array, axis=1)

    # Read storage_curve:
    storage_curve = np.array(pd.read_csv(sc_path, sep='\t'))
    time_seconds = np.array(time_to_seconds(df_time))

    # Interpolate initial values
    outflow_array[0, :] = calculate_outflow_rk(initial_wl)

    # Runge-Kutta method
    h_t = initial_wl
    for t in range(0, df_time.shape[0]-1):
        # k1 --------------------------------------------------------------------------------------------------------
        inflow_1 = interpolate_value(time_seconds[t], time_seconds, total_inflow)  # inflow t1

        # outflow as a function of water level (h_t): create function to calculate outflow, either by an equation (e.g.
        # for an orifice) or interpolating from a outflow vs wl curve (.txt file) EXAMPLE:
        outflow_1 = np.sum(calculate_outflow_rk(h_t))  # Outflow for the t1 time step

        # To calculate the area, calculate the volume for 'h', and then divide by 'h' (approximation)
        a_1 = interpolate_value(h_t, storage_curve[:, 0], storage_curve[:, 1])/h_t
        k1 = (inflow_1 - outflow_1)/a_1

        # k2 ---------------------------------------------------------------------------------------------------
        # Repeat procedure from before
        in_2 = interpolate_value(time_seconds[t] + dt/2, time_seconds, total_inflow)
        out_2 = np.sum(calculate_outflow_rk(h_t + (k1*dt/2)))
        a_2 = interpolate_value(h_t + (k1*dt/2), storage_curve[:, 0], storage_curve[:, 1])/(h_t + (k1*dt/2))
        k2 = (in_2 - out_2)/a_2

        # k3 ---------------------------------------------------------------------------------------------------
        # Repeat procedure from before
        in_3 = interpolate_value(time_seconds[t] + dt / 2, time_seconds, total_inflow)
        out_3 = np.sum(calculate_outflow_rk(h_t + (k2*dt / 2)))
        a_3 = interpolate_value(h_t + (k2 * dt/2), storage_curve[:, 0], storage_curve[:, 1]) / (h_t + (k2*dt/2))
        k3 = (in_3 - out_3) / a_3

        # k4 ---------------------------------------------------------------------------------------------------
        # Repeat procedure from before
        in_4 = interpolate_value(time_seconds[t] + dt, time_seconds, total_inflow)
        out_4 = np.sum(calculate_outflow_rk(h_t + (k3*dt)))
        a_4 = interpolate_value(h_t + (k3 * dt), storage_curve[:, 0], storage_curve[:, 1]) / (h_t + (k3*dt))
        k4 = k3*(in_4 - out_4) / a_4

        # h(t+1)
        h_t_1 = h_t + (k1 + 2*k2 + 2*k3 + k4)*dt / 6
        outflow_array[t+1, :] = calculate_outflow_rk(h_t_1)

        # For next iteration, if new h is higher than max water level, the excess water is removed, and the water level
        # is the max water level
        if h_t_1 > h_max:
            h_t_1 = h_max
        h_array[t+1, 0] = h_t_1

    # Total flow format:
    total_flow = np.c_[inflow_array, outflow_array]
    return total_flow


def monthly_inflow_avg(df_time, flow_array):
    """
    Function calculates the monthly inflow volume, by averaging the total inflow (for each sub-catchment individually)
    for a given month to get average monthly inflow (m3/s) for each, and then multiplying it by 3600 * 24 * days in
    the corresponding month to get total volume (m3) for each sub-catchment.

    Args:
    ----------------------------------------
    :param df_time: data frame DatetimeIndex or time Series, with original time intervals
    :param flow_array: np.array, with inflow data

    :return: data frame, with new time intervals, and np.array, with the monthly inflow volume for each month in the
    analysis time range.
    """
    # total_inflow = np.sum(flow_array, axis=1)
    month_dates, month_flow_total = resample_time(df_time, flow_array, interval=2)

    days_in_month = np.array(month_dates.days_in_month).reshape(month_flow_total.shape[0], 1)

    month_volume = np.multiply(month_flow_total, days_in_month) * 3600 * 24

    return month_dates, month_volume


# Sediment yield/concentration functions
def read_sediment_data(catchment_list):
    """
    Functions loops through each .txt file in the input sediment data folder and extracts the total sediment yield
    (ton/month) for each sub-catchment, as well as the date data. It then copies the SY data to an array, which must be
    in the following sub-catchment order: [Devoll, Holta, Zalli and Skebices], or in the order dictated by the variable
    'catchment_order' in config.py.

    Args:
    ----------------------------------
    :param catchment_list: list of strings, with the path to each .txt file in the input sediment folder

    :return: data frame, with monthly time intervals, corresponding to each SY data, and np.array, with the monthly
    total sediment yield for each sub-catchment (in each column), and for each time interval (row)

    Note: the name of the catchment, as it appears in the variable 'catchment_order', must be in the .txt file
    corresponding to that sub-catchment.

    """
    for i in range(0, len(catchment_list)):
        df = pd.read_csv(catchment_list[i], sep='\t')
        if i == 0:
            # Generate array to save data
            sy_array = np.full((df.shape[0], 4), 0.0)
            # Get time data
            date_df = pd.to_datetime(df['Date'], format='%Y%m')

        # read total_sy column
        temp_array = np.array(df['Total Sediment Yield [ton/month]'])
        # Assign SY data to corresponding column
        for c, c_name in enumerate(catchment_order):
            if c_name.lower() in catchment_list[i].lower():
                sy_array[:, c] = temp_array
                break
            else:
                if c == len(catchment_order)-1:
                    print(f"{catchment_list[i]} has no corresponding catchment in 'catchment_order' variable.")
    return date_df, sy_array


def calculate_concentration(sy_array, monthly_vol_array, sy_dates, vol_dates):
    """
    Function calculates the monthly volume concentration for each month in the discharge data, by following these steps:
    1. Trim the SY data to match the months considered in the discharge data
    2. Divide the SY data (ton/month) by the total monthly volume (m3/month), on a sub-catchment basis (column-wise)
    and then multiply by 1000 kg/ton to get the mass concentration (kg/m3)
    3. Divide the mass concentration (kg/m3) by the sediment density (kg/m3) to get the monthly volume concentration
    (m3/m3)

    Args:
    --------------------------------------------------------------------
    :param sy_array: np.array, with SY data (ton/month) for each sub-catchment
    :param monthly_vol_array: np.array, with the monthly volume data for each sub-catchment (m3/month)
    :param sy_dates: series, with dates considered in the SY data (each row is a different monthly interval)
    :param vol_dates: series, with dates considered in the discharge data, and which will be used in the model

    :return: np.array, with monthly volume concentration for each sub-catchment, for the time range in the discharge
    analysis, AND df, datetime index with trimmed dates (months) included in simulation.
    """
    # Filter sediment yield data for the months in the discharge data:
    df_total = pd.DataFrame(data=sy_array, index=sy_dates)
    df_trim = df_total.loc[vol_dates[0]:vol_dates[-1]]

    # Extract trimmed Sy data to array (ton/month)
    trimmed_sy_array = np.array(df_trim)
    trimmed_monthly_dates = df_trim.index

    # Calculate monthly mass concentration (kg/m3):
    mass_concent = np.divide(trimmed_sy_array, monthly_vol_array) * 1000
    # Calculate volume concentration (m3/m3)
    vol_concent = np.divide(mass_concent, sediment_density)

    return vol_concent, trimmed_monthly_dates


def build_concentration_timei(total_concentration, df_month, df_simulation_time):
    """
    Function build the timei file part with concentration values. There must be 4 entries per sub-catchment (with a
    total of 4 subcatchments): 1st is always 0 and the other 3 are 1/3 of the total monthly volume concentration.
    It first calculates 1/3 of the monthly volume concentration for each sub-catchment, and it then fills the concent.
    for each time interval in the simulation, by looping through each month in the simulation time, and then checking
    the month in each simulation time step (in order) and assigning the corresponding concentrations.

    Args:
    -----------------------------------------------
    :param total_concentration: np.array, with monthly volume concentrations for each sub-catchment (tx4 size)
    :param df_month: df, time series with each month for the volume concentration data (tx1)
    :param df_simulation_time: df, time series with each time interval in the final simulation (ix1, i>t)

    :return: np.array (concent_grain_fractions) with concentration values for each grain size fraction, and each
    sub-catchment for the final timei file (size ix16).
    """
    # Divide each volume concentration equally into 3 grain sizes:
    gs_concentration = total_concentration / 3

    concent_grain_fractions = np.full((df_simulation_time.shape[0], 4*4), 0.0)
    start_value = 0
    for m, month in enumerate(df_month):  # Loop through each month
        for i in range(start_value, df_simulation_time.shape[0]):  # Loop through each simulation time step
            interval = df_simulation_time[i]
            if month.month == interval.month and month.year == interval.year:
                pass
            else:
                # When month ends, fill all time steps in that month in the results array
                concent_grain_fractions[start_value:i, 1:4] = gs_concentration[m, 0]   # Devoll
                concent_grain_fractions[start_value:i, 5:8] = gs_concentration[m, 1]   # Holta
                concent_grain_fractions[start_value:i, 9:12] = gs_concentration[m, 2]  # Zalli
                concent_grain_fractions[start_value:i, 13:] = gs_concentration[m, 3]   # Skebices

                start_value = i  # To start in in next month, and avoid looping through all time steps
                break
    return concent_grain_fractions


# TIMEI FILES
def upstream_downstream_data(up, down, df_time):
    """
    Function fills the initial part (upstream and downstream data) of the 'I' dataset part of the timei file. It keeps
    the upstream and downstream discharges at 0, and assigns the upstream and downstream water level for each
    simulation time step.
    The order of columns is [upstream_discharge, downstream_discharge, upstream_wl, downstream_wl] or [0, 0, num, num]

    Args
    ---------------------------------
    :param up: float, with constant upstream water level (str, with path where .txt file with upstream WL data is)
    :param down: float, with constant downstream water level (str, with path where .txt file with downstream WL data is)
    :param df_time: df, time series with time steps in simulation

    :return: np.array, with (ix4) size array, with upstream and downstream discharges (value of 0) and water levels
    (constants).

    Note: it currently only works for a constant upstream and downstream water level. If a time-dependent water level
    is to be used, the corresponding code must be written to read the .txt file with the data, the time discretization,
    and writing it to the corresponding array.
    """
    us_ds_array = np.full((df_time.shape[0], 4), 0.0)
    if type(up) == int or float and type(down) == int or float:  # If input data are constant numerical type
        us_ds_array[:, 2] = upstream_wl
        us_ds_array[:, 3] = downstream_wl
    else:
        print("The program is not yet fit to read upstream and downstream water levels from .txt file")

    return us_ds_array


def build_timei_file(up_down_array, concent_array, flow_array, df_time):
    """
    Function generates the final timei file and saves it to a data frame variable. It then saves the timei file as a
    tab delimited text file (with no extension), called 'timei'.
    - Comprised of the "I" and "Q" data sets.
    - The "I" data set has the following format: [I, Time (s), 0, 0, upstream water level(m), downstream water
    level (m), conc. Devoll grain fraction 1 (m³/m³), conc. Devoll grain fraction 2, conc. Devoll grain fraction 3,
    conc. Devoll grain fraction 4 …] (ix22)
    - The "Q" data set has the following format: [Q Time(s)	Devoll [m³/s]	Holta [m³/s]	Zalli [m³/s]
    Skebices [m³/s]	Turbine inlet [m³/s]    Spillway [m³/s] 0   0   0]

    Args:
    --------------------------
    :param up_down_array: np.array, with upstream and downstream data (ix4) - from 'upstream_downstream_data' function
    :param concent_array: np.array, with volume concentration per grain size fractions and subcatchment (ix16) - from
    'build_concentration_timei' function
    :param flow_array: np.array, with inflow/outflow data (ix6) - from 'timei_total_flows' function
    :param df_time: df, time series in seconds, for each time step in final simulation

    :return: saves df (timei_df) with final boundary condition data for timei file
    """
    seconds_array = np.array(df_time).astype(int)

    # Letter df:
    i_letter_df = pd.DataFrame(data=np.full((seconds_array.shape[0], 1), "I"), columns=['Letter'])
    q_letter_df = pd.DataFrame(data=np.full((seconds_array.shape[0], 1), "Q"), columns=['Letter'])

    # time df
    time_df = pd.DataFrame(data=seconds_array, columns=['time'])

    # Combine upstream/downstrean data with concentrations for I data
    i_concat_array = np.concatenate((up_down_array, concent_array), axis=1)
    i_concat_df = pd.DataFrame(data=i_concat_array, columns=np.arange(1, i_concat_array.shape[1]+1))

    # Add 3 zeros at the end of flow data for Q data
    q_concat_array = np.concatenate((flow_array, np.full((seconds_array.shape[0], 3), 0.0)), axis=1)
    q_concat_df = pd.DataFrame(data=q_concat_array, columns=np.arange(1, q_concat_array.shape[1]+1))

    # Final I and Q df
    i_df = pd.concat([i_letter_df, time_df, i_concat_df], axis=1)
    q_df = pd.concat([q_letter_df, time_df, q_concat_df], axis=1)

    # TIMEI df format
    timei_df = pd.concat([i_df, q_df], axis=0)

    # Save timei
    timei_df.to_csv(os.path.join(resuls_folder, "timei"), header=False, index=False, sep='\t', na_rep="")


# Runge-Kutta method
def interpolate_value(value, x, y):
    """
    Function interpolates between two values from a np array. Currently, only a linear interpolation is done, but future
    modifications could include higher-order interpolations.

    Args:
    ----------------------------
    :param value: x value to interpolate for
    :param x: np.array (1D) with x values
    :param y: np.array (1D) with y values

    :return: interpolated y value, corresponding to the input x value
    """
    # Linear interpolation
    int_value = np.interp(value, x, y)
    return int_value


def calculate_outflow_rk(h):
    """
    Example of how to calculate outflow based on an outflow origice, using the orifice equation
    [q=pi/4*Cd*d^2*sqrt(2*g*h)], considering a diameter of 2m and an orifice coefficient Cd of 0.611. Fow flows above
    the one calculated for 173 m, the spillway begins to work such that the inflow equals the outflow, and the intake is
    located at a height of 161, such that below said height, no water leaves through the intake.

    Args
    ------------
    :param h: water level

    :return: np.array, with two values: the turbine outflow and the spillway outflow

    NOTE: This is just an example and thus must be modified to suit the needs of the project, e.g. using a
    stage-discharge curve or another equation, with different (appropriate) parameters

    NOTE: Instead of a maximum turbine capacity, one could set the maximum water level, which determines the max
    outflow through the turbine, and if h > h_max, then the extra height (h-h_max) is used to calculate the water
    through the spillway, for example. The h_max can be an input to the function, or just read from config directly
    """
    d = 1.15  # Diameter
    if h < 161:
        q_turbine = 0
        q_spillway = 0
    elif 161 <= h <= h_max:
        q_turbine = (math.pi/4)*0.611*math.pow(d, 2)*math.sqrt(2*9.81*h)
        q_spillway = 0
    else:
        q_turbine = (math.pi/4)*0.611*math.pow(d, 2)*math.sqrt(2*9.81*174)

        h_above_spillway = h - h_max
        q_spillway = 0  # Use equation for spillway, using h above spillway (Poleni equation)

    return np.array([q_turbine, q_spillway])

