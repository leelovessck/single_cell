#!/usr/bin/env python
# coding: utf-8

# In[2]:


# !unzip sc_PE_allcells_with_metadata_29-May-2023.txt.zip


# In[5]:


from tqdm import tqdm
file_path = "sc_PE_allcells_with_metadata_29-May-2023.txt"
line_list=[]
with open(file_path, "r") as file:
    for i, line in tqdm(enumerate(file)):
        if i<3:
            line_flag=line.split('\t')
            line_list.append(line_flag)
        else:
            break


# In[6]:


list1=['48-1', '53-1', '61-1', '63-1', '63-2', '67-1', '67-2', '68-1',
       '68-2', '77-1', '78-1', '86-1', '87-1', '88-1', '92-1', '99-1',
       '101-1', '103-1', '104-1', '105-1', '106-1', '107-1', '107-2',
       '108-1', '113-2', '114-2', '116-1', '117-1', '119-1', '120-1',
       '120-2', '120-2']


# In[7]:


for ii,sample in enumerate(list1):
    globals()['flag'+str(ii+1)]=[True if value==sample or i==0 else False for i,value in enumerate(line_flag)]


# In[8]:


out_path_1="sample_48-1.txt"
out_path_2="sample_53-1.txt"
out_path_3="sample_61-1.txt"
out_path_4="sample_63-1.txt"
out_path_5="sample_63-2.txt"
out_path_6="sample_67-1.txt"
out_path_7="sample_67-2.txt"
out_path_8="sample_68-1.txt"
out_path_9="sample_68-2.txt"
out_path_10="sample_77-1.txt"
out_path_11="sample_78-1.txt"
out_path_12="sample_86-1.txt"
out_path_13="sample_87-1.txt"
out_path_14="sample_88-1.txt"
out_path_15="sample_92-1.txt"
out_path_16="sample_99-1.txt"
out_path_17="sample_101-1.txt"
out_path_18="sample_103-1.txt"
out_path_19="sample_104-1.txt"
out_path_20="sample_105-1.txt"
out_path_21="sample_106-1.txt"
out_path_22="sample_107-1.txt"
out_path_23="sample_107-2.txt"
out_path_24="sample_108-1.txt"
out_path_25="sample_113-2.txt"
out_path_26="sample_114-2.txt"
out_path_27="sample_116-1.txt"
out_path_28="sample_117-1.txt"
out_path_29="sample_119-1.txt"
out_path_30="sample_120-1.txt"
out_path_31="sample_120-2.txt"


# In[9]:


with open(file_path, "r") as file,\
     open(out_path_1, "w") as output_file_1,\
     open(out_path_2, "w") as output_file_2,\
     open(out_path_3, "w") as output_file_3,\
     open(out_path_4, "w") as output_file_4,\
     open(out_path_5, "w") as output_file_5,\
     open(out_path_6, "w") as output_file_6,\
     open(out_path_7, "w") as output_file_7,\
     open(out_path_8, "w") as output_file_8,\
     open(out_path_9, "w") as output_file_9,\
     open(out_path_10, "w") as output_file_10,\
     open(out_path_11, "w") as output_file_11,\
     open(out_path_12, "w") as output_file_12,\
     open(out_path_13, "w") as output_file_13,\
     open(out_path_14, "w") as output_file_14,\
     open(out_path_15, "w") as output_file_15,\
     open(out_path_16, "w") as output_file_16,\
     open(out_path_17, "w") as output_file_17:
    file_list=[output_file_1, output_file_2, output_file_3, output_file_4, output_file_5, output_file_6, output_file_7, output_file_8, output_file_9, output_file_10, output_file_11, output_file_12, output_file_13, output_file_14, output_file_15, output_file_16, output_file_17]
    for i, line in tqdm(enumerate(file)):
        if line:
            line_split=line.split('\t')
            for jj in range(1,18):
                line1=[i for i,j in zip(line_split,globals()['flag'+str(jj)]) if j==True]
                file_list[jj-1].write("\t".join(line1) + "\n")          
        else:
            break


# In[17]:


with open(file_path, "r") as file,\
     open(out_path_18, "w") as output_file_18,\
     open(out_path_19, "w") as output_file_19,\
     open(out_path_20, "w") as output_file_20,\
     open(out_path_21, "w") as output_file_21,\
     open(out_path_22, "w") as output_file_22,\
     open(out_path_23, "w") as output_file_23,\
     open(out_path_24, "w") as output_file_24,\
     open(out_path_25, "w") as output_file_25,\
     open(out_path_26, "w") as output_file_26,\
     open(out_path_27, "w") as output_file_27,\
     open(out_path_28, "w") as output_file_28,\
     open(out_path_29, "w") as output_file_29,\
     open(out_path_30, "w") as output_file_30:
    file_list=[output_file_18, output_file_19, output_file_20, output_file_21, output_file_22, output_file_23, output_file_24, output_file_25, output_file_26, output_file_27, output_file_28, output_file_29, output_file_30]
    for i, line in tqdm(enumerate(file)):
        if line:
            line_split=line.split('\t')
            for jj in range(18,31):
                line1=[i for i,j in zip(line_split,globals()['flag'+str(jj)]) if j==True]
                file_list[jj-18].write("\t".join(line1))
                file_list[jj-18].write("\n")
        else:
            break

