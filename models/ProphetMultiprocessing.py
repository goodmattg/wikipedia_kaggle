
# coding: utf-8

# In[3]:

import pandas as pd
import numpy as np
from collections import defaultdict
from fbprophet import Prophet
import datetime
import multiprocessing
import sys


# # Read in Raw Data and Keys

# In[16]:

# Read in an take the transpose
rawData_df = pd.read_csv('../input/train.csv')
rawData_df.reset_index(inplace=True)
# Fill all NaN with zeros
rawData_df.fillna(value=0.0, inplace=True)
rawData_df.rename(columns={'index': 'mapId'}, inplace=True)


# In[5]:

keys_df = pd.read_csv('../input/key_1.csv')


# # Create Date Range

# In[7]:

date_start = datetime.date(2017, 1, 1)
date_end = datetime.date(2017, 3, 1)
num_days = 60

date_list = [date_start + datetime.timedelta(days=x) for x in range(0, num_days)]


# # Build a Mapping Dictionary: {Page -> Hash Id's}

# In[25]:

all_pages = rawData_df['Page'].tolist()

date_dict = dict(zip(keys_df['Page'].tolist(), keys_df['Id'].tolist()))


# In[9]:

mapping_dict = defaultdict(list)
for page in all_pages:
    for date in date_list:
        tmp = page + '_' + date.strftime('%Y-%m-%d')
        mapping_dict[page].append(date_dict[tmp])


# # Date Dataframe for Future Prediction

# In[43]:

future_df = pd.DataFrame(date_list, columns=['ds'])


# ### Set up Multiprocessing

# In[47]:

# Take the transpose
rawData_df = rawData_df.T


# In[45]:

FORECAST_DIR = 'forecasts/'


# In[51]:

def ProcessTimeSeries(idx):
    # Get the page name
    page_name = rawData_df.iloc[1, idx]

    try:

        test_df = rawData_df.iloc[2:, idx].to_frame().reset_index().fillna(method='ffill')
        test_df.columns = ['ds','y']
        test_df['ds'] = pd.to_datetime(test_df['ds'],format='%Y-%m-%d')

        # Train Out-of-the-box Prophet on the test dataframe
        m = Prophet(yearly_seasonality=True)
        m.fit(test_df)

        forecast = m.predict(future_df)
        forecast = forecast[['ds', 'yhat']]
        forecast['hash'] = pd.Series(mapping_dict[page_name], index=forecast.index)

        return {'page':page_name,
                'forecast':forecast.round(4)
        }

    except:
        with open(FORECAST_DIR + 'error_log', 'a') as f:
            f.write(page_name + '\n')
        f.closed
        return { 'page': page_name,
                'forecast':pd.DataFrame()
        }


# ### Run Multiprocessing

# In[52]:

pool = multiprocessing.Pool()

for resframe in pool.imap_unordered(ProcessTimeSeries, list(range(0,rawData_df.shape[1]))):
    with open(FORECAST_DIR + 'all_rows.csv', 'a') as f:
        resframe['forecast'].to_csv(f, header=False)

pool.close()
pool.join()



# In[ ]:



