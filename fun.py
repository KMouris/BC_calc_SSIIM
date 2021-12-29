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
    :param df_time: data frame, with original time intervals, in datetime format
    :param flow_array: np.array, with inflow data
    :param interval: 1, if re-sample is to daily, 2 if re-sample is to monthly

    :return: data frame, with new time intervals, and np.array, with the flow data averaged for the new time interval.

    Note: Function currently only re-samples to daily and monthly data, and averages the flow data for the new interval.
    If the user wants a different frequency, include options for the wanted frequency
    (check: https://regenerativetoday.com/pandas-data_range-function/)
    """

    df_total = pd.DataFrame(data=flow_array, index=df_time)

    if interval == 1:   # Day
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
    it to seconds.

    Args:
    ------------------------------------
    :param df_time: df, time series with each entry being a time interval for the simulation

    :return: df, time series in seconds.
    """
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


def monthly_inflow_avg(df_time, flow_array):
    """
    Function calculates the monthly inflow volume, by averaging the total inflow (for each sub-catchment individually)
    for a given month to get average monthly inflow (m3/s) for each, and then multiplying it by 3600 * 24 * days in
    the corresponding month to get total volume (m3) for each sub-catchment.

    Args:
    ----------------------------------------
    :param df_time: data frame, with original time intervals, in datetime format
    :param flow_array: np.array, with inflow data

    :return: data frame, with new time intervals, and np.array, with the monthly inflow volume for each month in the
    analysis time range.
    """
    # total_inflow = np.sum(flow_array, axis=1)
    month_dates, month_flow_total = resample_time(df_time, flow_array, interval=2)

    days_in_month = np.array(month_dates.days_in_month).reshape(month_flow_total.shape[0], 1)

    month_volume = np.multiply(month_flow_total, days_in_month) * 3600 * 24

    return month_dates, month_volume


# Soil yield/concentration functions
def read_soil_data(catchment_list):
    """
    Functions loops through each .txt file in the input soil data folder and extracts the total soil yield (ton/month)
    for each sub-catchment, as well as the date data. It then copies the SY data to an array, which must be in the
    following sub-catchment order: [Devoll, Holta, Zalli and Skebices], or in the order dictated by the variable
    'catchment_order' in config.py.

    Args:
    ----------------------------------
    :param catchment_list: list of strings, with the path to each .txt file in the input soil folder
    :return: data frame, with monthly time intervals, corresponding to each SY data, and np.array, with the monthly
    total soil yield for each subcatchment (in each column), and for each time interval (row)

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
    3. Divide the mass concentration (kg/m3) by the soil density (kg/m3) to get the monthly volume concentration (m3/m3)

    Args:
    --------------------------------------------------------------------
    :param sy_array: np.array, with SY data (ton/month) for each sub-catchment
    :param monthly_vol_array: np.array, with the monthly volume data for each sub-catchment (m3/month)
    :param sy_dates: series, with dates considered in the SY data (each row is a different monthly interval)
    :param vol_dates: series, with dates considered in the discharge data, and which will be used in the model

    :return: np.array, with monthly volume concentration for each sub-catchment, for the time range in the discharge
    analysis, AND df, datetime index with trimmed dates (months) included in simulation.
    """
    # Filter soil yield data for the months in the discharge data:
    df_total = pd.DataFrame(data=sy_array, index=sy_dates)
    df_trim = df_total.loc[vol_dates[0]:vol_dates[-1]]

    # Extract trimmed Sy data to array (ton/month)
    trimmed_sy_array = np.array(df_trim)
    trimmed_monthly_dates = df_trim.index

    # Calculate monthly mass concentration (kg/m3):
    mass_concent = np.divide(trimmed_sy_array, monthly_vol_array) * 1000
    # Calculate volume concentration (m3/m3)
    vol_concent = np.divide(mass_concent, soil_density)

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

