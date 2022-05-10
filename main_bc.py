"""
Program generates the boundary conditions for the SSIIMM2 numerical software, based on discharge data from WaSim
software and sediment yield calculations based on the following program:

Mouris, K., Morales Oreamuno, M.F., Schwindt, S., 2021a. Sediment_Load_Calculation.
https://github.com/KMouris/Sediment_Load_Calculation.

Edit: Surac, J. 2022
    Added seasonal changing water level with a threshold based water control logic.
    This changes the mass equation from inflow = outflow to storage = inflow - outflow.

"""

from fun import *

# Create results folder
if not os.path.exists(results_folder):
    os.makedirs(results_folder)

# READ INPUT FILES
# Read SY data:
filenames_sy = glob.glob(sy_folder + "/*.txt")

# read .b16 text file with discharge data to df

q_df = pd.read_csv(q_path, sep='\t', header=0, skiprows=[1, 2])
# rename to convert to date time
q_df.rename(columns={'YY': 'year', 'MM': 'month', 'DD': 'day', 'HH': 'hour'}, inplace=True)
# Convert first 4 columns to datetime and save to array
time_df = pd.to_datetime(q_df[['year', 'month', 'day', 'hour']])

# GET FLOW DATA
# Extract discharge data
q_array = extract_discharge(q_df)

# Option to convert to daily or monthly frequency
if time_interval != 0:
    time_df, q_array = resample_time(time_df, q_array, time_interval)

# Calculate outflows and water levels with a dynamic threshold based water control logic
timei_total_flows, timei_us_ds_array = calculate_outflows_seasonal_wl(q_array, time_df)

# GET CONCENTRATION DATA
# Calculate monthly volume from inflow data (for concentration data)
month_time_df, monthly_volume = monthly_inflow_avg(time_df, q_array)

# Read total sediment yield data and corresponding dates
sy_dates_df, sy_array = read_sediment_data(filenames_sy)

# Compare flow and concentration available dates, and trim to smallest date range
if sy_dates_df[0] > month_time_df[0] or sy_dates_df[sy_dates_df.shape[0] - 1] < month_time_df[-1]:
    print("Flow dates must be changed")
    timei_total_flows, time_df, monthly_volume, month_time_df, timei_us_ds_array = compare_flow_sediment_dates(
        sy_dates_df, time_df, timei_total_flows, monthly_volume, month_time_df, timei_us_ds_array)

# Calculate monthly concentration, for date range of discharge data:
total_concentration_array, trimmed_sy_dates = calculate_concentration(sy_array, monthly_volume, sy_dates_df,
                                                                      month_time_df)

# Separate monthly concentration for the 3 inflow grain size fractions (equally) for each sub-catchment
timei_concent_array = build_concentration_timei(total_concentration_array, trimmed_sy_dates, time_df)

# Build TIMEI file
timei_df = build_timei_file(timei_us_ds_array, timei_concent_array, timei_total_flows, time_df)
