import csv
import os

from selenium import webdriver
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait

chrome_option = webdriver.ChromeOptions()
chrome_option.add_argument(argument='--headless')
browser = webdriver.Chrome(options=chrome_option)   # options=chrome_option
wait = WebDriverWait(browser, 2)  # 最长等待时间为3S
# wait_v = WebDriverWait(browser, 20)
browser.maximize_window()


def crawler(g_name):
    if g_name == "na":
        return "unknown_type"

    base_url = "https://www.genecards.org/cgi-bin/carddisp.pl?gene="
    url = base_url + g_name
    browser.get(url)
    print("Searching:", g_name, end="...")
    try:
        g_type = wait.until(EC.presence_of_element_located(
            (By.CSS_SELECTOR, "#aliases_descriptions > div.row "
                              "> div.col-xs-12.col-md-9 > div:nth-child(2) > div > div"))).text
    except TimeoutException:
        try:
            temp = wait.until(EC.presence_of_element_located(
                (By.CSS_SELECTOR, "#cardPage > div:nth-child(1) > div > div > div.gc-section-header.gc-first "
                                  "> div.col-xs-8.col-md-10.col-sm-9.gene-symbol-description.row-no-padding "
                                  "> div > span.gc-category"))).text
            if temp == "Protein Coding":
                g_type = "protein_coding"
            elif temp == "Pseudogene":
                g_type = "pseudogene"
            else:
                print("----------------------------------", temp)
                g_type = temp

        except TimeoutException:
            g_type = "unknown_type"
    print("\t\ttype: ", g_type)
    return g_type


def combine_bed():
    all_bed = {}

    bed_path = "../results/circexplore2_output/"
    for _, dirs, _ in os.walk(bed_path):
        for dir_name in dirs:
            temp_path = bed_path + dir_name
            for _, _, filenames in os.walk(temp_path):
                for f in filenames:
                    file_path = temp_path + "/" + f
                    bed_file = open(file_path)
                    for line in bed_file.readlines():
                        info = line.strip().split()
                        rna_id = info[3]
                        rna_type = info[12]
                        rna_name = info[13]
                        if rna_id not in all_bed:
                            all_bed[rna_id] = [rna_type, rna_name]
    return all_bed


def load_circrna(all_bed):
    all_rna = {}
    diff_types = []
    rna_path = "../results/differential_expression/circRNA/"
    for _, dirs, _ in os.walk(rna_path):
        for dir_name in dirs:
            temp_path = rna_path + dir_name
            diff_types.append(dir_name)
            for _, _, filenames in os.walk(temp_path):
                for f in filenames:
                    if "whole_differential_expression.txt" in f:
                        file_path = temp_path + "/" + f
                        rna_file = open(file_path)
                        for line in rna_file.readlines()[1:]:
                            info = line.strip().split()
                            rna_id = info[0]
                            log2_fold_change = float(info[2])
                            p_value = float(info[5])
                            if rna_id not in all_rna:
                                all_rna[rna_id] = [log2_fold_change, p_value]
                            else:
                                all_rna[rna_id].extend([log2_fold_change, p_value])
                        break

    for each in all_rna.keys():
        if each in all_bed:
            all_rna[each].append(all_bed[each])
        # print(all_rna[each])
    return all_rna, diff_types


def load_gene_type():
    gene_type_file = open("./my_genbe_type.txt")
    gene_type = {}
    for line in gene_type_file.readlines():
        g_name, g_type = line.strip().split(" ")
        gene_type[g_name] = g_type
    gene_type_file.close()

    return gene_type


def generate_csv(all_rna, diff_types):
    combine_result = open("../results/combined_differential_expression.csv", "w")
    csv_title = "circrna_ID,"
    for each in diff_types:
        csv_title += ("log2FoldChange_%s,p_value_%s," % (each, each))
    csv_title += "circType,geneName,geneName2,geneName3,geneType\n"
    combine_result.write(csv_title)

    gene_type = load_gene_type()
    for each in all_rna.keys():
        value = all_rna[each]

        g_type = "unknown_type"
        if "," in value[6][1]:
            g_names = value[6][1].split(",")
            for g_name in g_names:
                g_name = g_name.lower()
                if g_name in gene_type:
                    g_type = gene_type[g_name]
                    break
                else:
                    g_type = crawler(g_name)
            if len(g_names) == 2:
                value[6][1] += ","
        else:
            g_name = value[6][1].lower()
            if g_name in gene_type:
                g_type = gene_type[g_name]
            else:
                g_type = crawler(g_name)
            value[6][1] += ",,"

        line = "%s,%lf,%lf,%f,%f,%f,%f,%s,%s,%s\n" \
               % (each, value[0], value[1], value[2], value[3], value[4], value[5], value[6][0], value[6][1], g_type)
        combine_result.write(line)

    combine_result.close()


if __name__ == '__main__':
    my_bed = combine_bed()
    my_rna, my_diff = load_circrna(my_bed)
    print("Circle RNA number: %d" % len(my_rna))

    generate_csv(my_rna, my_diff)
