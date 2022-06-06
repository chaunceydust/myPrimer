import re
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import os
from bs4 import BeautifulSoup
import numpy as np
import requests

# df3 = pd.DataFrame(columns=['Gene_id', 'Species', 'blast_strain_num', 'blast_target_num','protein_anno','primer_f', 'primer_r', 'primer_sp_num','primer_target_num'])
def blast(seq):
    gene_id = seq.split('\n')[0].strip('>')
    target = gene_id.split('_')[0].split(' ')[1]
    if not os.path.exists('result/'+target):
        os.mkdir('result/'+target)
    if gene_id in ''.join(os.listdir('result/'+target)):
        return
    org = 'bacteria (taxid:2)'
    df = pd.DataFrame(columns=['Accession', 'Title', 'Identity', 'Alignment length', 'Seq. start', 'Seq. stop', 'Gene'])
    s = Service(r"D:/2345Downloads/chromedriver.exe")
    browser = webdriver.Chrome(service=s)

    browser.get('https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?LINK_LOC=BlastHome')
    browser.find_element(By.XPATH, '//*[@id="PRIMER_PRODUCT_MIN"]').clear()
    browser.find_element(By.XPATH, '//*[@id="PRIMER_PRODUCT_MIN"]').send_keys(100)
    browser.find_element(By.XPATH, '//*[@id="ORGANISM"]').clear()
    browser.find_element(By.XPATH, '//*[@id="ORGANISM"]').send_keys(org)
    browser.find_element(By.XPATH, '//*[@id="PRIMER_SPECIFICITY_DATABASE"]/option[4]').click()
    browser.find_element(By.XPATH, '//*[@id="seq"]').send_keys(seq)
    browser.find_element(By.XPATH, '//*[@id="searchForm"]/div[3]/div[1]/input').click()
    # 显式等待
    wait = WebDriverWait(browser, 600)
    wait.until(EC.url_contains('job_key'))
    #wait.until(EC.staleness_of(browser.find_element(By.ID, 'statInfo')))
    table = browser.find_elements(By.XPATH, '//*[@id="descr"]/tbody/tr')
    count = 0
    for i in range(1, len(table)):
        line = table[i].find_elements(By.XPATH, 'td')
        df = df.append({'Accession': line[0].text, 'Title': line[1].text, 'Identity': line[2].text,
                        'Alignment length': line[3].text, 'Seq. start': line[4].text, 'Seq. stop': line[5].text,
                        'Gene': line[6].text}, ignore_index=True)
        if target not in line[1].text:
            count = count + 1
    sum = len(df)

    df.to_csv('result/%s/%s_%s_%s.csv' % (target, gene_id, sum, count), header=True, index=None)
    browser.find_element(By.XPATH, '//*[@id="userGuidedForm"]/div/div[1]/input').click()
    try:
        if 'No primers were found' in browser.find_element(By.XPATH, '//*[@id="content"]/div/ul/li/p').text:
            return
    except:
        pass
    wait = WebDriverWait(browser, 1200)
    wait.until(EC.url_contains('job_key'))
    url = browser.current_url
    with open('url.txt', "a") as file:
        file.write(gene_id+'\t'+target+'\t'+url+ "\n")

def primer(i):
    Ori_seq = i.split('\t')[0]
    target = i.split('\t')[1]
    url = i.split('\t')[2]
    if os.path.exists("result/%s/Primer_%s.txt" % (target, Ori_seq)):
        return
    r = requests.get(url)
    soup = BeautifulSoup(r.content, "lxml")

    data_list = []
    for i in soup.find('div', {'id': 'alignments'}).children:
        for n in i.find_all_next(attrs={"class": "prPairInfo"}):
            data_list.append(n.get_text())

    result_list = []
    for data in data_list:
        file = [x for x in data.split("\n") if x]
        Primer = file[0]
        FPrimer = [re.findall("(?<=Forward primer)(\w+)Plus", x)[0] for x in file if
                   re.findall("(?<=Forward primer)(\w+)Plus", x)][0]
        RPrimer = [re.findall("(?<=Reverse primer)(\w+)Minus", x)[0] for x in file if
                   re.findall("(?<=Reverse primer)(\w+)Minus", x)][0]
        all_species = [x.split(",")[0] for x in file if x.startswith(">")]
        target_species = [x for x in all_species if re.findall(target, x)]
        notarget_species = [x for x in all_species if not x in target_species]
        length = [int(re.findall("(?<=product length = )(\d+)", x)[0]) for x in file if
                  re.findall("(?<=product length = )(\d+)", x)]
        mean_len = round(np.mean(length), 3)
        std_length = round(np.std(length), 3)
        result_list.append(
            [Primer, Ori_seq, target, FPrimer, RPrimer, len(all_species), len(target_species), mean_len, std_length,
              ",".join(notarget_species)])

        df_out = pd.DataFrame(result_list,
                              columns=["Primer_Name", "Ori", "Target", "FP", "RP", "SP_Num", "Target_Num", "Mean_Len",
                                       "Std_Len", "Other_SP"]).drop_duplicates()
        df_out.to_csv("result/%s/Primer_%s.txt" % (target, Ori_seq), sep="\t", index=False)


with open('exp_5.fa','r')as f:
    all_fa = f.read().split('>')[1:]
all_faseq = list(map(lambda a:'>'+a,all_fa))
# seq='>gene_id_4367211\nATGGATTATAAAAAAATACTAGATCAGATCGTGAGTGGCGAATTAAAAGAATATAGAGTAAGTCCTAAAGATGCTTTTCAATTTCAATTAGCACTTAGAAGCTATGGTAAAAGGCAAAGTATAACTGGAATTGCTCAGCGTGGTGGGGATATAATTTATACCTTATCTGATAATGAAAATTAA'
# gene_id=seq.split('\n')[0].strip('>')
# target=gene_id.split('_')[0].split(' ')[1]

#map(primer,all_faseq)
for i in all_faseq:
    blast(i)

with open('url.txt','r')as f:
    all = f.read().split('\n')

for i in all:
    print(i)
    primer(i)

