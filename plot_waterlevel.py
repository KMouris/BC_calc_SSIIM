
from fun import *

import matplotlib.pyplot as mp
import time
from datetime import datetime


def plot_debug(res):
    df_plot = pd.read_csv(
        os.path.join(res, "debug.csv"),
        parse_dates=[0],
        header=None,
        sep='\t',
        names=["timestamp", "Water level", "Inflow", "Turbine", "Outflow"],
        converters={
            "timestamp": np.float64,
            "Water level": np.float64,
            "Inflow": np.float64,
            "Turbine": np.float64,
            "Outflow": np.float64
        }
    )

    df_plot["Date"] = [datetime.fromtimestamp(int(float(i))) for i in df_plot["timestamp"]]
    # df_plot = pd.DataFrame(data=plot_data, columns=)

    #df_plot.plot(x="Date", y="Water level", kind="line", style='-', figsize=(20, 9))
    df_plot.plot(x="Date", y=["Outflow", "Turbine"], kind="area", figsize=(20, 9))

    mp.show()


def plot_timei (res, filename):
    df_plot = pd.read_csv(
        os.path.join(res, filename),
        header=None,
        sep='\t',
        names=["state","timestamp","Water level"],
        usecols=[0,1,4],
        dtype={
            "state": str,
            "timestamp": np.float64,
            "Water level": np.float64
        }
    )
    # Only use input data
    df_plot = df_plot[df_plot["state"]=="I"]

    try:
        # Use start date given
        df_plot["Date"] = [datetime.fromtimestamp(int(i) + int((time.mktime(timei_date_start.timetuple())))) for i in df_plot["timestamp"]]
    except:
        # Use standard UNIX time to start (1970 01 01)
        df_plot["Date"] = [datetime.fromtimestamp(int(i)) for i in
                           df_plot["timestamp"]]
    #df_plot = pd.DataFrame(data=plot_data, columns=)

    #df_plot=df_plot.set_index("Date")
    ax = df_plot.plot(x="Date",y="Water level", kind="line", figsize=(50, 25))

    ax.set_ylabel('Wasserstand',fontdict={'fontsize':60})
    ax.set_xlabel('Datum',fontdict={'fontsize':60})
    ax.set_xticklabels(ax.get_xticks(),fontdict={'fontsize':50})
    ax.set_yticklabels(ax.get_yticks(),fontdict={'fontsize':50})

    mp.show()

res = results_folder
filename = "timei"
#filename = "timei_code_alt"

plot_timei(res, filename)
plot_debug(res)