#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd


# In[2]:


import matplotlib.pyplot as plot


# In[3]:


import numpy as np


# In[8]:


from sklearn.linear_model import LinearRegression


# In[9]:


mouse_info_df=pd.read_csv('Mouse_metadata.csv')
result_df=pd.read_csv('Study_results.csv')


# In[10]:


list_mouse_id=list(result_df.groupby(['Mouse ID','Timepoint'], as_index=False).filter(lambda y: len(y) > 1)['Mouse ID'].unique())


# In[11]:


###Summary Statistics Table###
for mouse_id in list_mouse_id:  
  result_df=result_df[result_df['Mouse ID']!=mouse_id]
print(result_df['Tumor Volume (mm3)'].describe() )


# In[12]:


full_df=pd.merge(mouse_info_df,result_df,on=['Mouse ID'])


# In[13]:



print(full_df[['Drug Regimen','Tumor Volume (mm3)']].groupby('Drug Regimen').describe())


# In[14]:


###Bar plot using datrame.plot()###
total_mice_per_regimen=mouse_info_df.groupby('Drug Regimen', as_index=False).size()


# In[15]:


drug_regimen_list=total_mice_per_regimen['Drug Regimen'].tolist()


# In[16]:


mice_count_list=total_mice_per_regimen['size'].tolist()


# In[17]:


bar_df = pd.DataFrame({'Drug Regimen':drug_regimen_list, 'Mice Count':mice_count_list})
bar_df.plot.bar(x='Drug Regimen', y='Mice Count', rot=70, title="Number of total mice for each treatment regimen")
plot.show()


# In[52]:


### Bar plot using matplotlib###
fig = plot.figure(figsize = (10, 5))
plot.bar(drug_regimen_list, mice_count_list, color ='blue',
        width = 0.4)
plot.xlabel("Drug Regimen")
plot.ylabel("Mice Count")
plot.title("Number of total mice for each treatment regimen")
plot.show()


# In[19]:





# In[20]:


### Pie plot using datrame.plot()####
gender_count=mouse_info_df.groupby('Sex', as_index=False).size()
print(gender_count)


# In[21]:


gender_list=gender_count['Sex'].tolist()


# In[22]:


gender_count_list=gender_count['size'].tolist()


# In[23]:


pie_df = pd.DataFrame({'Sex':gender_list, 'Sex Count':gender_count_list}, index=['Female','Male'])
pie_df.plot.pie(y='Sex Count', autopct='%1.1f%%', figsize=(5, 5), title="Distribution of Female or Male mice")
plot.show()


# In[53]:


### Pie plot using matplotlib####
fig = plot.figure(figsize = (10, 5))
plot.pie(gender_count_list,autopct='%1.1f%%',labels =['Female','Male'])
plot.legend(title = "Sex:")
plot.title("Distribution of Female or Male mice")
plot.show()


# In[51]:





# In[28]:


###using datrame.plot()####
tumor_volume_list=result_df[result_df['Mouse ID']=='x401']['Tumor Volume (mm3)'].tolist()


# In[29]:


timepoint_list=result_df[result_df['Mouse ID']=='x401']['Timepoint'].tolist()


# In[30]:


line_df = pd.DataFrame({'Tumor Volume (mm3)':tumor_volume_list, 'Timepoint':timepoint_list})
line_df.plot.line(x='Tumor Volume (mm3)', y='Timepoint', title="Tumour Volume vs. Time Point for Mouse ID=x401 treated with Capomulin")
plot.show()


# In[32]:


###using datrame.plot()####
tumor_avg_per_weight=full_df[full_df['Drug Regimen']=='Capomulin'][['Weight (g)','Tumor Volume (mm3)']].groupby('Weight (g)', as_index=False).mean()
 


# In[33]:


mouse_weight_list=tumor_avg_per_weight['Weight (g)'].tolist()


# In[34]:


average_tumour_volume_list=tumor_avg_per_weight['Tumor Volume (mm3)'].tolist()


# In[35]:


scatter_df = pd.DataFrame({'Mouse Weight':mouse_weight_list, 'Average Tumour Volume':average_tumour_volume_list})
scatter_df.plot.scatter(x='Mouse Weight', y='Average Tumour Volume',c='DarkBlue', title="Mouse Weight versus Average Tumour Volume for the Capomulin treatment regimen")
plot.show()


# In[36]:


pearsoncorr=scatter_df.corr(method='pearson')
print("correlation coefficient mouse weight and average tumour volume for the Capomulin treatment:")
print(pearsoncorr)


# In[37]:


###Linear regression####


# In[38]:


linear_regressor = LinearRegression()
X=scatter_df['Mouse Weight'].values.reshape(-1,1)
Y=scatter_df['Average Tumour Volume'].values.reshape(-1,1)
linear_regressor.fit(X,Y )
#predicting the set result
Y_pred = linear_regressor.predict(X)
print(Y_pred)
plot.scatter(X, Y, color = 'blue')
plot.plot(X, Y_pred, color = 'red')
plot.title('Mouse Weight versus Average Tumour Volume for the Capomulin treatment regimen')
plot.xlabel('Mouse Weight')
plot.ylabel('Average Tumour Volume')
plot.show()


# In[39]:


most_promising_regimens=['Capomulin', 'Ramicane', 'Infubinol', 'Ceftamin']
max_timepoint=full_df[full_df['Drug Regimen'].isin(most_promising_regimens)].groupby(['Mouse ID'], as_index=False).Timepoint.max()


# In[40]:


final_tumor_vol_per_mouse=pd.DataFrame(columns = ['Mouse ID', 'Drug Regimen', 'Final Tumor Volume (mm3)'])


# In[41]:


for index,rows in max_timepoint.iterrows():
  get_df=full_df[(full_df['Mouse ID']==rows['Mouse ID']) & (full_df['Timepoint']==rows['Timepoint'])]  
  final_tumor_vol_per_mouse=final_tumor_vol_per_mouse.append({'Mouse ID' : get_df['Mouse ID'].values[0],'Drug Regimen':get_df['Drug Regimen'].values[0], 'Final Tumor Volume (mm3)':get_df['Tumor Volume (mm3)'].values[0]}, ignore_index = True)


# In[42]:


print(final_tumor_vol_per_mouse) 


# In[43]:


data=[]


# In[44]:


stats=pd.DataFrame(columns = ['Drug Regimen','IQR','First quartile','Third quartile'])
for i in most_promising_regimens:  
  data.append(final_tumor_vol_per_mouse[final_tumor_vol_per_mouse['Drug Regimen']==i]['Final Tumor Volume (mm3)'] )
  finaltumor_vol_per_regimen=final_tumor_vol_per_mouse[final_tumor_vol_per_mouse['Drug Regimen']==i]['Final Tumor Volume (mm3)']
  q25 = finaltumor_vol_per_regimen.quantile(0.25)
  q75 = finaltumor_vol_per_regimen.quantile(0.75) 
  iqr = q75 - q25
  stats=stats.append({'Drug Regimen':i,'IQR':iqr,'First quartile':q25, 'Third quartile':q75},ignore_index = True)


# In[45]:


print(stats)


# In[46]:


####box plot###


# In[47]:


fig = plot.figure(figsize = (10, 5))
plot.boxplot(data, flierprops=dict(markerfacecolor='g', marker='D'))
plot.xticks([1, 2, 3, 4], most_promising_regimens,rotation=10)
plot.title("Final tumour volume for all four treatment regimens")
plot.show()


# In[ ]:




