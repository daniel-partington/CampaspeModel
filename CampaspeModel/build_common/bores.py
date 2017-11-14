import matplotlib.pyplot as plt

def plot_bores_large_stdev_errors(bores_obs_time_series):
    # Looking at bores with large standard deviations for errors:
    b = bores_obs_time_series[bores_obs_time_series['value'] != 0.0]
    b = b[b['active'] == True]
    c = b[['name', 'value']].groupby('name').std()
    c.hist()
    d = c[c['value'] > 4]
    for bore in d.index.tolist():
        b_df = bores_obs_time_series[bores_obs_time_series['name'] == bore]
        b_df_min = b_df['datetime'].min()
        b_df_max = b_df['datetime'].max()
        b_df_sd = 4.0 * b_df['value'].std()
        b_df_mean = b_df['value'].mean()
        ax = b_df.plot(x='datetime', y='value', marker='o', label=bore)
        plt.plot([b_df_min, b_df_max], [b_df_mean, b_df_mean], label='mean')
        plt.plot([b_df_min, b_df_max], [b_df_mean + b_df_sd, b_df_mean + b_df_sd], label='2.5$\sigma$')
        plt.plot([b_df_min, b_df_max], [b_df_mean - b_df_sd, b_df_mean - b_df_sd], label='2.5$\sigma$') 