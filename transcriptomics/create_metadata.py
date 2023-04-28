import pandas as pd 
import re 
import os

dir_name = "PRJNA867309"

path = "/Users/davide/Desktop/transcriptomics/" + dir_name

files_dir = path +"/files/"

with open(files_dir+'biosample_result.txt') as f:
    file = f.read()
    print(file)

# extract features from the biosample file and put them into a dictionary

dic = {}
features = list(set([match.group().strip('/=') for match in re.finditer(r'/.*=', file)]))
for feature in features:
    dic[feature] = [re.search('".*"', match.group()).group().strip('"') for match in re.finditer(f'/{feature}.*', file)]

print(dic)

# In[3]: create a dataframe with the features

df = pd.DataFrame.from_dict(dic)


# In[4]:


runinfo = pd.read_csv(files_dir + "SraRunInfo.csv")


# In[5]:


run = runinfo["Run"]


# In[7]: insert if the "Run" column and Bioproject


df.insert(0, "Run", run)

for num in range(len(run)):
	bioproject.append(dir_name)

df.insert(0, "Bioproject", bioproject)



# In[8]:


df.to_csv(files_dir + "metadata.tsv", sep ="\t", index = False)
