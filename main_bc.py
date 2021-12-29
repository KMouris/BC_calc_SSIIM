"""
Program generates the boundary conditions for the SSIIMM2 numerical software, based on discharge data from WaSim
software and soil yield calculations based on the following program:

Mouris, K., Morales Oreamuno, M.F., Schwindt, S., 2021a. Sediment_Load_Calculation.
https://github.com/KMouris/Sediment_Load_Calculation.
"""

from fun import *

# Create results folder
if not os.path.exists(resuls_folder):
    os.makedirs(resuls_folder)

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
timei_total_flows = calculate_outflows_constant_wl(q_array, turbine_capacity)

# ADD 3 ZEROS AT THE END OF THE TOTAL FLOWS (always) AND A 'Q' AT THE BEGINNING IN FINAL DF

# GET CONCENTRATION DATA ----------------------------------------------------------------------------------------
# Calculate monthly volume form inflow data(for concentration data)
month_time_df, monthly_volume = monthly_inflow_avg(time_df, q_array)

# Read total soil yield data and corresponding dates
sy_dates_df, sy_array = read_soil_data(filenames_soil)

# Calculate monthly concentration, for date range of discharge data:
total_concentration_array, trimmed_sy_dates = calculate_concentration(sy_array, monthly_volume, sy_dates_df,
                                                                      month_time_df)

# Separate monthly concentration for the 3 inflow grain size fractions (equally) for each sub-catchment
timei_concent_array = build_concentration_timei(total_concentration_array, trimmed_sy_dates, time_df)

# upstream and downstream WL data for timei file
timei_us_ds_array = upstream_downstream_data(upstream_wl, downstream_wl, time_df)

# Time steps in seconds
seconds_df = time_to_seconds(time_df)

# Build TIMEI file
build_timei_file(timei_us_ds_array, timei_concent_array, timei_total_flows, seconds_df)

x=1

