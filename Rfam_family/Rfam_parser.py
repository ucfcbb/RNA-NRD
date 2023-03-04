#!/usr/bin/env python
# coding: utf-8

# In[224]:


from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By

from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.support.ui import Select
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException

import time 


# In[225]:


#DRIVER_PATH = 'D:\SOFTWARE\chromedriver_win32\chromedriver.exe' #For my laptop
DRIVER_PATH = '/media/nabila/NABILA6TB/chromedriver_linux64/chromedriver' #For my work pc

#DRIVER_PATH = Service('/media/nabila/NABILA6TB/chromedriver_linux64/chromedriver')
# In[226]:


#driver = webdriver.Chrome(executable_path=DRIVER_PATH)
#driver = webdriver.Chrome(service=DRIVER_PATH)
driver.maximize_window() # For maximizing window
driver.implicitly_wait(20) # gives an implicit wait for 20 seconds



# In[227]:


options = Options()
options.headless = True
driver.get("https://rfam.org/search?q=entry_type:%22Family%22%20AND%20has_3d_structure:%22Yes%22")
SCROLL_PAUSE_TIME = 0.5


# In[228]:



# Family_name = open("Family_name.txt", "w")
#xpath = "/html/body/div[5]/div[4]/div[1]/div[2]/div[1]/div[2]/div/div/div/div[2]/div[2]/div[2]/button[contains(@class,'btn btn-default load-more col-md-3')]"
#xpath = "/html/body/div[5]/div[4]/div[1]/div[2]/div[1]/div[2]/div/div/div/div[2]/div[2]/div[2]/button"
xpath = "/html/body/div[4]/div[4]/div[1]/div[2]/div[1]/div[2]/div/div/div/div[2]/div[2]/div[2]/button" ### "Lod more" button
page_no = 1

### To find absolute xpath, click the element and then go to inspect/copy/full xpath

Family_list = {}
flag = 0

while(True):
    #lists = driver.find_element(By.XPATH, "/html/body/div[4]/div[4]/div[1]/div[2]/div[1]/div[2]/div/div/div/div[2]/div[2]/div[2]/ul")
    #print(lists)
    #rows =driver.find_element(By.XPATH, "/html/body/div[4]/div[4]/div[1]/div[2]/div[1]/div[2]/div/div/div/div[2]/div[2]/div[2]/ul/li[1]/div/div/h4/a") ### Each family name
    #print(rows)

    #/html/body/div[4]/div[4]/div[1]/div[2]/div[1]/div[2]/div/div/div/div[2]/div[2]/div[2]/ul/li[1]/div/div/h4/a

    lists = driver.find_element_by_xpath("/html/body/div[4]/div[4]/div[1]/div[2]/div[1]/div[2]/div/div/div/div[2]/div[2]/div[2]/ul") ### All the family boxes
    rows =lists.find_elements_by_xpath(".//li/div/div/h4/a") ### Each family name
    
    line_no = 1
    
   
    if(flag==1):
        line_no = 1
        for row in rows:
            print(line_no, ':\t', row.text)
#             Family_name.write(str(line_no) + ':\t' + str(row.text) + '\n')
            Family_list[line_no] = row.text
            line_no += 1        
        break    
         
   
    for row in rows:
        current_row = row.text
        
#         print(line_no, '    :::    ', row.text)
#         Family_name.write(str(line_no) + '    :::    ' + str(row.text) + '\n')
#         line_no += 1
    
    cls = driver.find_element_by_xpath(xpath).get_attribute("class")

    if cls == "btn btn-default load-more col-md-3":
        driver.execute_script("arguments[0].click();", WebDriverWait(driver, 20).until(EC.element_to_be_clickable((By.CSS_SELECTOR, "button[class='btn btn-default load-more col-md-3']"))))
        SCROLL_PAUSE_TIME = 0.5
    else:
        line_no = 1
        for row in rows:
            current_row = row.text
            
#             print(line_no, ':\t', row.text)
#             Family_name.write(str(line_no) + ':\t' + str(row.text) + '\n')
#             Family_list[line_no] = row.text
#             line_no += 1        
        # break
        flag = 1
        continue
    
    
# Family_name.close()


# In[229]:



#xpath = "/html/body/div[5]/div[4]/div[1]/div[2]/div[1]/div[2]/div/div/div/div[2]/div[2]/div[2]/button[contains(@class,'btn btn-default load-more col-md-3')]"
#xpath = "/html/body/div[5]/div[4]/div[1]/div[2]/div[1]/div[2]/div/div/div/div[2]/div[2]/div[2]/button"
xpath = "/html/body/div[4]/div[4]/div[1]/div[2]/div[1]/div[2]/div/div/div/div[2]/div[2]/div[2]/button"

while(True):
    #lists = driver.find_element_by_xpath("/html/body/div[5]/div[4]/div[1]/div[2]/div[1]/div[2]/div/div/div/div[2]/div[2]/div[2]/ul")
    lists = driver.find_element_by_xpath("/html/body/div[4]/div[4]/div[1]/div[2]/div[1]/div[2]/div/div/div/div[2]/div[2]/div[2]/ul")
    rows =lists.find_elements_by_xpath(".//li/div/div/div/div[2]/div[2]/div/ul/li[3]/span/a")
    line_no = 1
    link_list = []
   
    for row in rows:
        link = row.get_attribute("href")
        link_list.append(link)
#         print(line_no, ': ', link)
        line_no += 1

    cls = driver.find_element_by_xpath(xpath).get_attribute("class")

    if cls == "btn btn-default load-more col-md-3":
        driver.execute_script("arguments[0].click();", WebDriverWait(driver, 20).until(EC.element_to_be_clickable((By.CSS_SELECTOR, "button[class='btn btn-default load-more col-md-3']"))))
    else:       
        break
        


# In[230]:


#pdb_family_list = open("PDB_family_name.txt", "w")
pdb_family_list = open("PDB_family_name.txt", "w")

for i in range(0, len(link_list)):

    cur_link = link_list[i]
    print("Line_no:", i+1, ':\t', cur_link)
    pdb_family_list.write("Line_no:" + str(i+1) + ':\t' + str(cur_link) + '\n')
    # driver.get("https://rfam.xfam.org/family/RF03013#tabview=tab6")
    driver.get(cur_link)
    driver.implicitly_wait(10)
    SCROLL_PAUSE_TIME = 0.5
    
    try:
        #table = driver.find_element_by_xpath("/html/body/div[5]/div[5]/div[1]/div[2]/div[7]/div[2]/div/table")
        table = driver.find_element_by_xpath("/html/body/div[4]/div[5]/div[1]/div[2]/div[7]/div[2]/div/table") ### Table inside "3D structures" page for each family
        rows =table.find_elements_by_xpath(".//tr")
        count = 0

        for row in rows:
            count += 1
            if(count == 1):
                continue

            rowTds = row.find_elements_by_xpath(".//td")
            column = []
            for tds in rowTds:
                column.append(tds.text)
                #print(tds.text)

            print(str(column[1]) + '\t' + str(column[2]) + '\t' + str(Family_list[i+1]))
            pdb_family_list.write(str(column[1]) + '\t' + str(column[2]) + '\t' + str(Family_list[i+1]) + '\n')
    except:
        continue

pdb_family_list.close()




# In[ ]:




