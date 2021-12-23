"""
Program generates the boundary conditions for the SSIIMM2 numerical software, based on discharge data from WaSim
software and soil yield calculations based on the following program:

Mouris, K., Morales Oreamuno, M.F., Schwindt, S., 2021a. Sediment_Load_Calculation.
https://github.com/KMouris/Sediment_Load_Calculation.
"""

from config import *


# FUNCTIONS ----------------------------------------------------------------------------------------------------
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


def calculate_outflows_constant_wl(inflow_array, turbine_capacity):
    """
    Function calculates the outflows based on a constant water level, meaning that inflows = outflows. The outflows are
    through a turbine, with a maximum capacity, and a spillway, in case the inflow exceeds the turbine capacity. The
    total inflow is calculated as the sum of the inflows from each sub-catchment (row-wise sum), and the outflow is
    calculated for each time interval (each row).

    Args:
    --------------------------------------------
    :param inflow_array: np.array, with the inflow data for each sub-catchment (in this case 4). Each column is a
    different sub-catchment and each row corresponds to a different time interval.
    :param turbine_capacity: float, turbine capacity

    :return: np.array with [inflows (4), turbine outflow, spillway outflow], in that order.
    """
    total_inflow = np.sum(inflow_array, axis=1)

    turbine_outflow = np.where(total_inflow > turbine_capacity, turbine_capacity, total_inflow)
    spillway_outflow = np.where(total_inflow > turbine_capacity, total_inflow - turbine_capacity, 0)

    # Total flow data
    total_flow = np.c_[inflow_array, turbine_outflow, spillway_outflow]

    return total_flow


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
    analysis.
    """
    # Filter soil yield data for the months in the discharge data:
    df_total = pd.DataFrame(data=sy_array, index=sy_dates)
    df_trim = df_total.loc[vol_dates[0]:vol_dates[-1]]

    # Extract trimmed Sy data to array (ton/month)
    trimmed_sy_array = np.array(df_trim)

    # Calculate monthly mass concentration (kg/m3):
    mass_concent = np.divide(trimmed_sy_array, monthly_vol_array) * 1000
    # Calculate volume concentration (m3/m3)
    vol_concent = np.divide(mass_concent, soil_density)

    return vol_concent


# READ INPUT FILES ---------------------------------------------------------------------------------------------
# Read SY data:
filenames_soil = glob.glob(sy_folder + "/*.txt")

# read .b16 file with discharge data to df
q_df = pd.read_csv(q_path, sep='\t', header=0, skiprows=[1, 2])
# rename to convert to date time
q_df.rename(columns={'YY': 'year', 'MM': 'month', 'DD': 'day', 'HH': 'hour'}, inplace=True)
# Convert first 4 columns to datetime and save to array
time_df = pd.to_datetime(q_df[['year', 'month', 'day', 'hour']])

# GET FLOW DATA ----------------------------------------------------------------------------------------------
# Extract discharge data
q_array = extract_discharge(q_df)

# Option to convert to daily or monthly frequency
if time_interval != 0:
    time_df, q_array = modify_time_interval(time_df, q_array, time_interval)
    resample_time(time_df, q_array, time_interval)

# Calculate outflows and save inflows and outflows to array
total_flows = calculate_outflows_constant_wl(q_array, turbine_capacity)

# ADD 3 ZEROS AT THE END OF THE TOTAL FLOWS (always) AND A 'Q' AT THE BEGINNING IN FINAL DF

# GET CONCENTRATION DATA ----------------------------------------------------------------------------------------
# Calculate monthly volume form inflow data(for concentration data)
month_time_df, monthly_volume = monthly_inflow_avg(time_df, q_array)

# Read total soil yield data and corresponding dates
sy_dates_df, sy_array = read_soil_data(filenames_soil)

# Calculate monthly concentration, for date range of discharge data:
total_concentration_array = calculate_concentration(sy_array, monthly_volume, sy_dates_df, month_time_df)

# SEPARATE MONTHLY CONCENTRATION FOR THE 3 INFLOW GRAIN SIZE FRACTIONS, EQUALLY

# CALCULATE TIME IN SECONDS, TO A LIST

# CONSTANTS FOR 'I' PART OF TIMEI
x=1

