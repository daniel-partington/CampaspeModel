import calendar

import pandas as pd


def year_period_2_datetime(year, period):
    year_lower = int(year[0].split("/")[0])
    year_upper = int(year[0].split("/")[1])
    if year_upper > 90:
        year_upper = 1900 + year_upper
    else:
        year_upper = 2000 + year_upper
    # end if
    first_fin = ['Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    date = "{} {}".format(year_lower if period[0] in first_fin else year_upper, period[0])
    return pd.to_datetime(date)


def get_diversions(fname, plot=False, custom_summary=False):
    '''
    Function to parse all the diversions data for Campaspe River
    from excel spreadsheet
    '''
    Campaspe_diversions = pd.read_excel(fname, sheetname=['Campaspe River Use', 'WWC', 'CID Diversions'])
    CD = Campaspe_diversions
    # Process Campaspe Private Diversions first (Ml/d)
    CD['Campaspe River Use'].columns = CD['Campaspe River Use'].iloc[0]
    CD['Campaspe River Use'] = CD['Campaspe River Use'].dropna(axis=1, how='all')
    CD['Campaspe River Use'] = CD['Campaspe River Use'].reindex(CD['Campaspe River Use'].index.drop(0))
    CD['Campaspe River Use'] = CD['Campaspe River Use'].drop('Cumulative Total (ML)', axis=1)
    CD['Campaspe River Use'] = CD['Campaspe River Use'].drop(CD['Campaspe River Use'].columns[[0, 1]], axis=1)
    CD['Campaspe River Use'] = CD['Campaspe River Use'].iloc[[0, 1, 2, 3, 6, 7, 10, 11]]
    CD['Campaspe River Use'].index = CD['Campaspe River Use'].apply(
        lambda row: year_period_2_datetime(row['Year'], row['Period']), axis=1)
    CD['Campaspe River Use'] = CD['Campaspe River Use'].drop(['Year', 'Period'],  axis=1)
    Reaches = [x for x in CD['Campaspe River Use'].iloc[0] if type(x) != float]
    CD['Campaspe River Use'] = CD['Campaspe River Use'].drop('Reach', axis=1)
    CD['Campaspe River Use'].columns = Reaches
    index_month2days = CD['Campaspe River Use'].apply(
        lambda row: calendar.monthrange(row.name.year, row.name.month)[1], axis=1)
    for col in CD['Campaspe River Use'].columns:
        CD['Campaspe River Use'][col] = CD['Campaspe River Use'][col] / index_month2days

    CD['Campaspe River Use'] = CD['Campaspe River Use'].astype(float)
    if plot:
        ax = CD['Campaspe River Use'].plot()
        ax.set_ylabel("Diversion [Ml/d]")

    # Process Waranga Western Channel next (Ml/d)
    CD['WWC'] = CD['WWC'].reindex(CD['WWC'].index.drop(0))
    CD['WWC'].index = CD['WWC']['Site Name']
    CD['WWC'] = CD['WWC'].drop('Site Name', axis=1)
    CD['WWC'] = CD['WWC'].astype(float)
    if plot:
        ax2 = CD['WWC'].plot()
        ax2.set_ylabel("Exchange [Ml/d]")

    # Process CID (Ml/d)
    CD['CID Diversions'].colummns = CD['CID Diversions'].iloc[0]
    CD['CID Diversions'].columns = [x for x in CD['CID Diversions'].iloc[0]]
    CD['CID Diversions'] = CD['CID Diversions'].reindex(CD['CID Diversions'].index.drop(0))
    CD['CID Diversions'].index = CD['CID Diversions']['Date']
    CD['CID Diversions'] = CD['CID Diversions'].drop('Date', axis=1)
    CD['CID Diversions'] = CD['CID Diversions'].astype(float)
    if plot:
        ax3 = CD['CID Diversions'].plot()
        ax3.set_ylabel("Diversion [Ml/d]")

    # Summarising of data:
    Total_diversions = CD['Campaspe River Use'].sum(axis=1)
    WWC_in = CD['WWC']['WWC TO CAMPASPE (RO317)'] + CD['WWC']['WWC TO CAMPASPE (RO317).1']
    WWC_out = CD['WWC']['CAMPASPE PUMP TO WWC -REG'] + CD['WWC']['CAMPASPE PUMP TO WWC -UNREG']
    CID_total = CD['CID Diversions'].sum(axis=1)
    if custom_summary:
        if plot:
            ax1 = Total_diversions.plot(label='Private diversions')
            CID_total.plot(ax=ax1, label='CID')
            WWC_in.plot(ax=ax1, label='WWC in')
            WWC_out.plot(ax=ax1, label='WWC out')
            ax1.set_title("Major diversions into and out of the Campaspe River")
            ax1.legend()
            ax1.set_xlabel("")
            ax1.set_ylabel("[Ml/d]")

    return CD


if __name__ == "__main__":
    # fname = r"C:\Workspace\part0075\MDB modelling\Campaspe_Data\SW\Campaspe_System_Data.xlsx"
    fname = r"C:\UserData\takuyai\ownCloud\campaspeim_data\Groundwater\2017-08-21\GW_data\Campaspe_data\SW\Campaspe_System_Data.xlsx"
    CD = get_diversions(fname, plot=True, custom_summary=True)
