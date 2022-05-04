"""
Module contains functions to read and calculate discharge (inflow/outflow) and volume concentrations for the generation
of the input boundary condition file for SSIIM2

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

    :param df_time: data frame, with original time intervals, in datetime format
    :param flow_array: np.array, with inflow data
    :param interval: 1, if re-sample is too daily, 2 if re-sample is to monthly

    :return: data frame, with new time intervals, and np.array, with the flow data averaged for the new time interval.

    Note: Function currently only re-samples to daily and monthly data, and averages the flow data for the new interval.
    If the user wants a different time interval or a different way to manage the flow data, they must do the
    corresponding changes.
    """

    start_date = df_time[0]
    end_date = df_time.iloc[-1]

    if interval == 1:  # Daily (frequency is 'D')
        freq_str = "D"
    elif interval == 2:  # monthly (frequency is "MS")
        freq_str = "MS"
    else:
        print(f"No resampling was done, since {interval} is not a valid entry")

        return df_time, flow_array

    resampled_dates = pd.to_datetime(pd.date_range(start_date, end_date, freq=freq_str).strftime("%Y%m%d").tolist(),
                                     format="%Y%m%d")
    resampled_flow = np.full((len(resampled_dates), flow_array.shape[1]), 0.0)
    initial = 0
    for d, date in enumerate(resampled_dates):  # Loop through each day to consider
        for i in range(initial, df_time.shape[0]):  # Loop through the original flow/date data
            if ((interval == 1) & (date.day != df_time[i].day)) | ((interval == 2) & (date.month != df_time[i].month)):
                resampled_flow[d, :] = np.mean(flow_array[initial:i, :], axis=0)

                initial = i
                break

    return resampled_dates, resampled_flow


def resample_time(df_time, flow_array, interval):
    """
    Function re-samples the initial time intervals to a new interval and averages the inflow data, for each watershed,
    for the new time interval.

    It first converts the flow array to a data frame, with the original time data as index. It then uses the df.resample
    tool to resample the data to a given frequency, by averaging the values.

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

    if interval == 1:  # Day
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

    :param df_time: DataFrame Series or DatetimeIndex, time series with each entry being a time interval for the
    simulation.

    :return: df, time series in seconds.

    Note: if df_time is a DatetimeIndex type, it is transformed to Series type.
    """
    if isinstance(df_time, pd.DatetimeIndex):  # if Datetime index, change to Series
        df_time = df_time.to_series()
    seconds = [np.int64((df_dt - df_time[:1].item()).total_seconds()) for df_dt in df_time]
    return seconds


# Discharge functions
def extract_discharge(input_df):
    """
    Function reads the raw .b16 discharge output file, in data frame format, and extracts the discharges for the Devoll,
    Holta, Zalli and Skebices sub-catchments. In the input file, columns with names {2, 3, 4, 5} (or indexes
    {6, 7, 8, 9}) corresponds to the catchments {Zalli, Skebices, Holta, Devoll}. The output file needs the following
    sub-catchment order: {Devoll, Holta, Zalli and Skebices}

    :param input_df:  data frame, with the original input data from the .b16 file

    :return: np array, with the discharges for the [Devoll, Holta, Zalli and Skebices] subcatchments, in that order
    """
    output_q = np.full((input_df.shape[0], 4), 0.0)

    output_q[:, 0] = input_df['5']  # Devoll
    output_q[:, 1] = input_df['4']  # Holta
    output_q[:, 2] = input_df['2']  # Zalli
    output_q[:, 3] = input_df['3']  # Skebices
    # output_q[:, 4] = input_df['0']  # Skebices

    return output_q


def calculate_outflows_seasonal_wl(inflow_array, df_time):
    """
    Function calculates the outflows and water levels based on a changing seasonal water level, meaning that
    inflows - outflows = storage. The total inflow is calculated as the sum of the inflows from each sub-catchment
    (row-wise sum), and the outflow and water level is calculated for each time interval (each row). The water level is
    linearly interpolated based on the current volume and a water level volume correlation file which is written in the
    variable q_storage.
    The water level is controlled with a dynamic threshold based logic. If the water level sinks 2 meters below the
    target level the storage mode is activated and the turbine has no output. If the water level is within (-2) - (-0)
    of the target level only 50 % of the inflow runs through the turbine. If the water level is within (+0) - (+0.1) of
    the target level the turbine runs with the inflow and only as high as the turbine capacity. When above the turbine
    runs always at maximum capacity. If the water level exceeds the 175 m water level an additional spillway drops the
    inflow.

    Note: This function is less accurate with time_interval = 1 (averaging by day) or time_interval = 2 (averaging by
    month). This is due to the different modis being applied to longer time frames. For example if the water level
    drops below the -2 a.s.l and the time_interval is set to two (monthly), all inflow within one month will be stored.

    :param inflow_array: np.array, with the inflow data for each sub-catchment (in this case 4). Each column is a
    different sub-catchment and each row corresponds to a different time interval.
    :param df_time: data frame DatetimeIndex or time Series, with original time intervals

    :return: np.array with [inflows (4), turbine outflow, spillway outflow], in that order.
    :return: np.array with [0, 0, downstream water level, upstream water level], in that order.
    """

    # Read file with wl and volume correlation
    q_vol_wl_correlation = pd.read_csv(q_storage, sep='\t', header=0)
    total_inflow = np.sum(inflow_array, axis=1)

    # Create empty array with target wl and volume
    us_ds_array = np.full((df_time.shape[0], 4), 0.0)
    plot_data = np.zeros(shape=(df_time.shape[0], 6))

    # Create empty array with outflow for every time step
    turbine_outflow = np.full((df_time.shape[0], 1), 0.0)
    spillway_outflow = np.full((df_time.shape[0], 1), 0.0)

    # Iterate through each data point
    for i, (df_dt, inflow) in enumerate(zip(df_time, total_inflow)):

        overflow = 0
        current_turbine_capacity = turbine_capacity

        # Apply target wl based on hydrological year
        target_wl = wl_threshold[df_dt.month - 1]

        # Fill upper and lower boundary for storage controll
        target_wl_upper = target_wl + target_wl_upper_boundary
        target_wl_lower_emergency = target_wl + target_wl_lower_boundary

        target_volume = q_vol_wl_correlation[q_vol_wl_correlation["Elevation"] == target_wl]["Volume"].values[0]

        # Apply current water level and volume
        if i == 0:
            current_volume = target_volume
            current_wl = target_wl
        else:
            # Waterlevel for current volumes
            lower_wl = q_vol_wl_correlation[q_vol_wl_correlation["Volume"] <= current_volume].values[-1]
            if lower_wl[0] == target_wl_maximum:
                upper_wl = lower_wl
            else:
                upper_wl = q_vol_wl_correlation[q_vol_wl_correlation["Volume"] >= current_volume].values[0]

            # Interpolate current water level
            if upper_wl[0] == lower_wl[0]:
                current_wl = upper_wl[0]
            else:
                current_wl = (current_volume - lower_wl[1]) / (upper_wl[1] - lower_wl[1]) * (
                        upper_wl[0] - lower_wl[0]) + lower_wl[0]

        # Water level is below emergency threshold therefore store all inflow
        if current_wl < target_wl_lower_emergency:
            current_turbine_capacity = 0

        # Water level is lower than the target water level therefore throttle turbine
        elif current_wl < target_wl:
            current_turbine_capacity = inflow * .5

        # Water level is above 175 m asl therefore spill additional inflow
        elif current_wl >= target_wl_maximum:
            if current_turbine_capacity < inflow:
                overflow = inflow - current_turbine_capacity

        # Water level is above the target threshold therefore use full turbine capacity
        elif current_wl > target_wl_upper:
            current_turbine_capacity = turbine_capacity

        # Water level is within threshold therefore keep steady
        else:
            if turbine_capacity >= inflow:
                current_turbine_capacity = inflow

        # New waterlevel according to mass equation
        # Inflow, turbine and overflow is [m³/s] while every time step is 3 hours. Therefore, multiply by 3600*3.

        if time_interval == 0:
            time_step = 3
        elif time_interval == 1:
            time_step = 24
        elif time_interval == 2:
            time_step = 24 * calendar.monthrange(df_dt.year, df_dt.month)[1]
        else:
            raise Exception(
                f"Error in calculate outflow function. Time interval is not valid: {time_interval}")

        current_volume = current_volume + (inflow - current_turbine_capacity - overflow) * 3600 * time_step

        # Output of the reservoir
        turbine_outflow[i] = current_turbine_capacity
        spillway_outflow[i] = overflow

        # Total water level data
        us_ds_array[i, 2] = current_wl
        us_ds_array[i, 3] = current_wl

        plot_data[i, 0] = np.int64((round(df_dt.timestamp())))
        plot_data[i, 1] = inflow
        plot_data[i, 2] = (inflow - current_turbine_capacity - overflow)
        plot_data[i, 3] = current_turbine_capacity
        plot_data[i, 4] = overflow
        plot_data[i, 5] = current_wl

    if plot_outflow_data:
        plot_inflow_outflow(pd.DataFrame(data=plot_data,
                                         columns=["timestamp", "inflow", "delta volume", "turbine", "overflow",
                                                  "water level"]))

    # Total flow data
    total_flow = np.c_[inflow_array, turbine_outflow, spillway_outflow]

    return total_flow, us_ds_array


def monthly_inflow_avg(df_time, flow_array):
    """
    Function calculates the monthly inflow volume, by averaging the total inflow (for each sub-catchment individually)
    for a given month to get average monthly inflow (m3/s) for each, and then multiplying it by 3600 * 24 * days in
    the corresponding month to get total volume (m3) for each sub-catchment.

    :param df_time: data frame DatetimeIndex or time Series, with original time intervals
    :param flow_array: np.array, with inflow data

    :return: data frame, with new time intervals, and np.array, with the monthly inflow volume for each month in the
    analysis time range.
    """
    month_dates, month_flow_total = resample_time(df_time, flow_array, interval=2)

    days_in_month = np.array(month_dates.days_in_month).reshape(month_flow_total.shape[0], 1)

    month_volume = np.multiply(month_flow_total, days_in_month) * 3600 * 24

    return month_dates, month_volume


def compare_flow_sediment_dates(sy_dates, df_time, timei_flow_array, monthly_vol_array, vol_dates, timei_us_ds_array):
    """
    Function compares the dates in the flow and the sediment data, and determines the date ranges included in both
    data sets, and selects the smallest date range, and trims the flow and monthly volume data arrays.

    Args:
    ------------------------------
    :param sy_dates: series, with dates considered in the SY data (each row is a different monthly interval)
    :param df_time: data frame DatetimeIndex or time Series, with original time intervals
    :param timei_flow_array: np.array, with inflow and outflow data, in timei format
    :param monthly_vol_array: np.array, with the monthly volume data for each sub-catchment (m3/month)
    :param vol_dates: series, with dates considered in the discharge data, and which will be used in the model
    :param timei_us_ds_array: np.array, with water levels upstream and downstream (m a.s.l)

    :return: modified flow and volume data (np.array) and dates (time Series or DatetimeIndex).

    Note: This function is called if the sediment date range is smaller than the flow date range. If the other way
    around, the sediment date range is modified in the "calculate_concentration" function.
    """
    # Determine start date:
    if sy_dates[0] > vol_dates[0]:
        # print("Start date is set by the sy data.")
        start_date = sy_dates[0]
    else:
        # print("Start date is set by the flow data.")
        start_date = vol_dates[0]

    # Determine end date:
    if sy_dates[sy_dates.shape[0] - 1] < vol_dates[-1]:
        # print("End date is set by the sy data.")
        end_date = sy_dates[sy_dates.shape[0] - 1]
    else:
        # print("End date is set by the flow data.")
        end_date = vol_dates[-1]

    # Get last day of month, last hour.
    end_date = end_date.replace(day=calendar.monthrange(end_date.year, end_date.month)[1])
    end_date = end_date.replace(hour=23, minute=59)

    if vol_dates[0] != start_date or vol_dates[-1] != end_date:
        trim_vol_array, trim_vol_date = trim_data_to_date(monthly_vol_array, vol_dates, start_date, end_date)
        trim_flow_array, trim_flow_date = trim_data_to_date(timei_flow_array, df_time, start_date, end_date)
        trim_timei_us_ds_array, trim_flow_date = trim_data_to_date(timei_us_ds_array, df_time, start_date, end_date)

        return trim_flow_array, trim_flow_date, trim_vol_array, trim_vol_date, trim_timei_us_ds_array
    else:
        return timei_flow_array, df_time, monthly_vol_array, vol_dates, timei_us_ds_array


# Sediment yield/concentration functions
def read_sediment_data(catchment_list):
    """
    Functions loops through each .txt file in the input sediment data folder and extracts the total sediment yield
    (ton/month) for each sub-catchment, as well as the date data. It then copies the SY data to an array, which must be
    in the following sub-catchment order: [Devoll, Holta, Zalli and Skebices], or in the order dictated by the variable
    'catchment_order' in config.py.

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
                if c == len(catchment_order) - 1:
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
    total of 4 subcatchments): 1st is always 0 and the other 3 are 1/3 of the total monthly volume concentration. It
    first calculates 1/3 of the monthly volume concentration for each sub-catchment, and it then fills the
    concentration for each time interval in the simulation, by looping through each month in the simulation time,
    and then checking the month in each simulation time step (in order) and assigning the corresponding concentrations.

    :param total_concentration: np.array, with monthly volume concentrations for each sub-catchment (tx4 size)
    :param df_month: df, time series with each month for the volume concentration data (tx1)
    :param df_simulation_time: df, time series with each time interval in the final simulation (ix1, i>t)

    :return: np.array (concent_grain_fractions) with concentration values for each grain size fraction, and each
    sub-catchment for the final timei file (size ix16).
    """
    # Divide each volume concentration equally into 3 grain sizes:
    gs_concentration = total_concentration / 3

    concent_grain_fractions = np.full((df_simulation_time.shape[0], 4 * 4), 0.0)
    start_value = 0
    for m, month in enumerate(df_month):  # Loop through each month
        for i in range(start_value, df_simulation_time.shape[0]):  # Loop through each simulation time step
            interval = df_simulation_time[i]
            if month.month == interval.month and month.year == interval.year:
                pass
            else:
                # When month ends, fill all time steps in that month in the results array
                concent_grain_fractions[start_value:i, 1:4] = gs_concentration[m, 0]  # Devoll
                concent_grain_fractions[start_value:i, 5:8] = gs_concentration[m, 1]  # Holta
                concent_grain_fractions[start_value:i, 9:12] = gs_concentration[m, 2]  # Zalli
                concent_grain_fractions[start_value:i, 13:] = gs_concentration[m, 3]  # Skebices

                start_value = i  # To start in in next month, and avoid looping through all time steps
                break
    # Fill in for last month
    concent_grain_fractions[start_value:, 1:4] = gs_concentration[m, 0]  # Devoll
    concent_grain_fractions[start_value:, 5:8] = gs_concentration[m, 1]  # Holta
    concent_grain_fractions[start_value:, 9:12] = gs_concentration[m, 2]  # Zalli
    concent_grain_fractions[start_value:, 13:] = gs_concentration[m, 3]  # Skebices

    return concent_grain_fractions


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

    :param up_down_array: np.array, with upstream and downstream data (ix4) - from 'upstream_downstream_data' function
    :param concent_array: np.array, with volume concentration per grain size fractions and subcatchment (ix16) - from
    'build_concentration_timei' function
    :param flow_array: np.array, with inflow/outflow data (ix6) - from 'timei_total_flows' function
    :param df_time: df, time series in seconds, for each time step in final simulation

    :return: saves df (timei_df) with final boundary condition data for timei file
    """

    mask_timeframe = get_timeframe_mask(df_time)

    df_time = df_time[mask_timeframe]

    # Apply time frame as a mask
    masked_df_time = time_to_seconds(df_time)
    masked_up_down_array = up_down_array[mask_timeframe]
    masked_concent_array = concent_array[mask_timeframe]
    masked_flow_array = flow_array[mask_timeframe]

    seconds_array = np.array(masked_df_time).astype(int)

    # Letter df:
    i_letter_df = pd.DataFrame(data=np.full((seconds_array.shape[0], 1), "I"), columns=['Letter'])
    q_letter_df = pd.DataFrame(data=np.full((seconds_array.shape[0], 1), "Q"), columns=['Letter'])

    # time df
    time_df = pd.DataFrame(data=seconds_array, columns=['time'])

    # Combine upstream/downstream data with concentrations for I data
    i_concat_array = np.concatenate((masked_up_down_array, masked_concent_array), axis=1)
    i_concat_df = pd.DataFrame(data=i_concat_array, columns=np.arange(1, i_concat_array.shape[1] + 1))

    # Add 3 zeros at the end of flow data for Q data
    q_concat_array = np.concatenate((masked_flow_array, np.full((seconds_array.shape[0], 3), 0.0)), axis=1)
    q_concat_df = pd.DataFrame(data=q_concat_array, columns=np.arange(1, q_concat_array.shape[1] + 1))

    # Final I and Q df
    i_df = pd.concat([i_letter_df, time_df, i_concat_df], axis=1)
    q_df = pd.concat([q_letter_df, time_df, q_concat_df], axis=1)

    # TIMEI df format
    timei_df = pd.concat([i_df, q_df], axis=0)

    # Save timei
    timei_df.to_csv(os.path.join(results_folder, "timei"), header=False, index=False, sep='\t', na_rep="")

    return timei_df


# General
def get_timeframe_mask(df_time):
    """
    Function creates a mask based on the given start date (dt_start) and end date (dt_end). This is only done if the
    restrict_timei_date tag is set to true.
    If there is no data found within the mask. Then the full mask is set to True and therefore no filter is
    applied.

    :param df_time: np.array, with the datetime object where the filter is applied to
    :return: pd.Series, with the shape of df_time and a boolean value for each entry.
    """
    true_mask = pd.Series([True] * df_time.shape[0])

    # Select specific time frame if stated
    if restrict_timei_date:
        dt_start = parser.parse(timei_date_start)
        dt_end = parser.parse(timei_date_end)

        mask_timeframe = ((df_time >= dt_start) & (df_time <= dt_end))

        # Check if there is valid data within the timeframe
        if mask_timeframe.sum() == 0:
            print(f"There was no data in the selected time frame. All available data was selected to process.")
            mask_timeframe = true_mask
        return mask_timeframe


def trim_data_to_date(data_array, date_df, start_date, end_date):
    """
    Function trims a data array and the corresponding dates to a different start and end date

    Args:
    -------------------------------
    :param data_array: np.array, with data to trim
    :param date_df: DataFrame, as time Series or DatetimeIndex
    :param start_date: Timestamp, with start date
    :param end_date: Timestamp, with end date

    :return: trimmed np.array (with data) and DatetimeIndex (with Timestamp dates)
    """
    df_total = pd.DataFrame(data=data_array, index=date_df)
    df_trim = df_total.loc[start_date:end_date]

    trimmed_array = np.array(df_trim)
    trimmed_dates = df_trim.index

    return trimmed_array, trimmed_dates


def plot_inflow_outflow(df_plot):
    """
    Function plots water level and the mass bilance in two plots. The dataframe has to have the following columns:
    "timestamp": int, with the unix time
    "water level": float, with the water levels in meter
    "inflow": float, with the total water inflow in m³/s
    "overflow": float, with the water dumped in the overflow in m³/s
    "turbine": float, with the water powering the turbine in m³/s

    The two plots are stored in the results_folder given in the config file as "waterlevel.png" and "massbilance.png"

    :param df_plot: pd.DataFrame, with the timestamp, water levels and inflow/outflow values
    """

    df_plot["Date"] = [datetime.datetime.fromtimestamp(int(float(i))) for i in df_plot["timestamp"]]
    # df_plot = pd.DataFrame(data=plot_data, columns=)

    # Trim data within mask
    mask = get_timeframe_mask(df_plot["Date"])

    df_plot = df_plot[mask]

    df_plot = df_plot.set_index("Date")

    fig1, ax1 = mp.subplots(figsize=[50, 25])
    fig2, ax2 = mp.subplots(figsize=[50, 25])

    # Water level plot
    df_plot.plot(ax=ax1, y="water level", label="Water level", kind="line", stacked=True, legend=False)

    ax1.set_title("Water level", fontsize=70)

    ax1.grid(b=True, which='major', axis='y', color='lightgrey', linestyle='-')
    ax1.grid(b=True, which='major', axis='x', color='lightgrey', linestyle='-')

    ax1.set_ylabel('Water level [m]', fontdict={'fontsize': 60})
    ax1.set_xlabel('Date ', fontdict={'fontsize': 60})
    ax1.tick_params(axis='both', labelsize=50)

    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels, loc='upper right', fontsize=60)

    fig1.savefig(os.path.join(results_folder, "waterlevel.png"))
    fig1.show()

    # Mass bilance plot
    df_plot.plot(ax=ax2, y="inflow", label="Inflow", color="blue", kind="line", linewidth=5, legend=False)
    df_plot.plot(ax=ax2, y=["turbine", "overflow"], label=["Turbine", "Overflow"], color=["green", "red"], kind="area",
                 stacked=True, legend=False)

    ax2.set_title("Mass bilance", fontsize=70)

    ax2.grid(b=True, which='major', axis='y', color='lightgrey', linestyle='-')
    ax2.grid(b=True, which='major', axis='x', color='lightgrey', linestyle='-')

    ax2.set_ylabel('Water throughput [m³/s]', fontdict={'fontsize': 60})
    ax2.set_xlabel('Date ', fontdict={'fontsize': 60})
    ax2.tick_params(axis='both', labelsize=50)

    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles, labels, loc='upper right', fontsize=60)

    fig2.savefig(os.path.join(results_folder, "massbilance.png"))
    fig2.show()
